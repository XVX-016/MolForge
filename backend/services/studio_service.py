
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

import httpx
import asyncio

class StudioService:
    def __init__(self):
        # Use global settings
        api_key = settings.GEMINI_KEY or settings.GOOGLE_API_KEY
        key_log = f"{api_key[:15]}...{api_key[-5:]}" if api_key else "NONE"
        print(f"DEBUG: StudioService key: {key_log}")
        
        if not api_key:
            logger.warning("GEMINI_KEY not found in environment variables. AI features will be disabled.")
            self.api_key = None
        else:
            self.api_key = "".join(api_key.split())
            logger.info("StudioService initialized with Gemini v1 REST API")

    def _extract_json(self, response_data: dict) -> dict:
        """Extract JSON from Gemini REST response"""
        try:
            text = response_data['candidates'][0]['content']['parts'][0]['text']
            # Find JSON block if model wrapped it in markdown
            start = text.find("{")
            end = text.rfind("}")
            if start != -1 and end != -1:
                return json.loads(text[start:end+1])
            return json.loads(text)
        except (KeyError, IndexError, json.JSONDecodeError) as e:
            logger.error(f"Failed to parse Gemini response: {e}")
            raise ValueError(f"Invalid response format from AI: {e}")

    def _is_safe_prompt(self, prompt: str) -> bool:
        """Relaxed safety filter"""
        return True

    async def process_command(self, prompt: str, context: Dict[str, Any], mode: str) -> Dict[str, Any]:
        """
        Processes a natural language command using Gemini v1 REST API.
        """
        if not self.api_key:
            return {
                "type": "NO_OP",
                "reason": "Gemini API key is missing. Please configure GEMINI_KEY in the backend environment."
            }

        try:
            # Prepare the request
            model_name = getattr(settings, "GEMINI_MODEL", "gemini-1.5-flash")
            url = f"https://generativelanguage.googleapis.com/v1/models/{model_name}:generateContent?key={self.api_key}"
            
            context_str = json.dumps(context)
            system_instruction = STUDIO_SYSTEM_PROMPT
            
            payload = {
                "contents": [
                    {
                        "role": "user",
                        "parts": [
                            {"text": f"System Instruction: {system_instruction}\n\nMode: {mode}\nContext: {context_str}\n\nUser Prompt: {prompt}"}
                        ]
                    }
                ],
                "generationConfig": {
                    "temperature": 0.0,  # Set to 0 for most predictable JSON
                    "topP": 0.95,
                    "topK": 40,
                    "maxOutputTokens": 1024
                }
            }

            async with httpx.AsyncClient(timeout=30.0) as client:
                response = await client.post(url, json=payload)
                
                if response.status_code != 200:
                    error_detail = response.text
                    logger.error(f"Gemini API Error ({response.status_code}): {error_detail}")
                    return {
                        "type": "NO_OP",
                        "reason": f"AI service error ({response.status_code}): {error_detail}"
                    }
                
                data = response.json()
                return self._extract_json(data)
                
        except Exception as e:
            logger.error(f"StudioService Error: {type(e).__name__}: {e}", exc_info=True)
            return {
                "type": "NO_OP",
                "reason": f"AI service processing error: {str(e)}"
            }

_studio_service: Optional[StudioService] = None

def get_studio_service() -> StudioService:
    global _studio_service
    if _studio_service is None:
        _studio_service = StudioService()
    return _studio_service
