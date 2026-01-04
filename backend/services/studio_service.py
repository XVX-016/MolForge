
import os
import logging
import json
from typing import Dict, Any, Optional
from pathlib import Path
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

# Try to import Gemini, but handle gracefully if not installed
try:
    import google.generativeai as genai
    GEMINI_AVAILABLE = True
except ImportError:
    GEMINI_AVAILABLE = False
    logger.warning("google-generativeai not installed. Studio service will not work.")

STUDIO_SYSTEM_PROMPT = """
You are MolForge Studio Architect, a high-precision molecular design orchestrator.
You ONLY output valid JSON. No conversational text.

Current Capabilities:
1. CREATE_MOLECULE: { "atoms": [], "bonds": [] }
2. ADD_ATOM: { "element": "C", "position": [x, y, z] }
3. REPLACE_ATOM: { "atomId": "id", "newElement": "N" }
4. REMOVE_ATOM: { "atomId": "id" }
5. ADD_BOND: { "from": "id1", "to": "id2", "order": 1 }
6. REMOVE_BOND: { "bondId": "id" }
7. OPTIMIZE_GEOMETRY: {}
8. SIMULATE_REACTION: {}
9. NO_OP: { "reason": "why" }

Rules:
- Default to NO_OP if a request is biologically/physically nonsensical.
- For CREATE_MOLECULE, use a sensible geometry if not specified.
- Use id format like 'a1', 'a2' for atoms and 'b1', 'b2' for bonds.
- Always include a "reason" field for scientific justification.
- If the user asks for property prediction or complex analysis, point them to the ANALYZE mode or provide a brief high-level summary in "reason".

Output Format:
{
  "type": "ACTION_TYPE",
  "payload": { ... },
  "reason": "..."
}
"""

from backend.config import settings

class StudioService:
    def __init__(self):
        if not GEMINI_AVAILABLE:
            raise ValueError("google-generativeai package not installed")

        # Use global settings
        api_key = settings.GEMINI_KEY or settings.GOOGLE_API_KEY
        key_log = f"{api_key[:15]}...{api_key[-5:]}" if api_key else "NONE"
        print(f"DEBUG: StudioService key: {key_log}")
        print(f"DEBUG: Key length: {len(api_key) if api_key else 0}")
        
        if not api_key:
            logger.warning("GEMINI_KEY not found in environment variables. AI features will be disabled.")
            self.api_key = None
            self.model = None
        else:
            self.api_key = "".join(api_key.split())
            logger.info(f"Gemini key loaded: {self.api_key[:6]}...")
            
            genai.configure(api_key=self.api_key)
            # Use simplest possible initialization that worked in script
            self.model = genai.GenerativeModel(model_name="models/gemini-1.5-flash")
            logger.info("StudioService initialized with models/gemini-1.5-flash (STABLE)")

    def _extract_json(self, text: str) -> dict:
        """Safely extract JSON from response text"""
        try:
            # First try direct parsing
            return json.loads(text)
        except json.JSONDecodeError:
            # Fallback to finding braces
            start = text.find("{")
            end = text.rfind("}")
            if start == -1 or end == -1:
                raise ValueError("No JSON object found in response")
            return json.loads(text[start:end+1])

    def _is_safe_prompt(self, prompt: str) -> bool:
        """
        Check if prompt is chemistry/studio related.
        Relaxed to allow more natural language interactions.
        """
        return True

    async def process_command(self, prompt: str, context: Dict[str, Any], mode: str) -> Dict[str, Any]:
        """
        Processes a natural language command into a structured action.
        """
        if not self.model:
            # Debug info for user
            env_path = Path(__file__).resolve().parent.parent / ".env"
            cwd = os.getcwd()
            has_key = "GEMINI_KEY" in os.environ or "GOOGLE_API_KEY" in os.environ
            
            return {
                "type": "NO_OP",
                "reason": f"Gemini API key is missing. Searched: {env_path} (Exists: {env_path.exists()}). CWD: {cwd}. EnvVars: {has_key}. Please configure GEMINI_KEY in the backend environment."
            }

        if not self._is_safe_prompt(prompt):
            return {
                "type": "NO_OP",
                "reason": "I can only assist with molecular design and chemistry-related tasks."
            }

        try:
            # Re-configure just in case
            genai.configure(api_key=self.api_key)
            
            # Using sync generate_content because async is currently returning 404 for this key
            model = genai.GenerativeModel(
                model_name="models/gemini-1.5-flash",
                system_instruction=STUDIO_SYSTEM_PROMPT
            )
            
            # Prepare context-aware prompt
            context_str = json.dumps(context)
            full_prompt = f"Mode: {mode}\nContext: {context_str}\n\nUser: {prompt}"
            
            import asyncio
            # Use asyncio.to_thread to run sync generation without blocking
            response = await asyncio.to_thread(
                model.generate_content,
                full_prompt,
                generation_config={
                    "temperature": 0.2,
                    "top_p": 0.95,
                    "top_k": 40,
                    "max_output_tokens": 1024,
                    "response_mime_type": "application/json"
                }
            )
            
            text = response.text
            return self._extract_json(text)
                
        except Exception as e:
            import traceback
            logger.error(f"StudioService Error Traceback: {traceback.format_exc()}")
            logger.error(f"StudioService Error: {type(e).__name__}: {e}")
            
            error_msg = str(e)
            if "401" in error_msg or "403" in error_msg or "API_KEY" in error_msg:
                user_msg = "AI service authentication failed. Please check API key."
            elif "429" in error_msg:
                user_msg = "Rate limit reached. Please try again later."
            else:
                user_msg = f"AI service error: {str(e)}"
            
            return {
                "type": "NO_OP",
                "reason": user_msg
            }

_studio_service: Optional[StudioService] = None

def get_studio_service() -> StudioService:
    global _studio_service
    if _studio_service is None:
        _studio_service = StudioService()
    return _studio_service
