
import os
import logging
import json
import re
from typing import Dict, Any, Optional
import httpx

logger = logging.getLogger(__name__)

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

class StudioService:
    def __init__(self):
        api_key = os.getenv("GEMINI_KEY")
        if not api_key:
            logger.warning("GEMINI_KEY not found in environment variables. AI features will be disabled.")
            self.api_key = None
        else:
            self.api_key = "".join(api_key.split())
            logger.info(f"Gemini key loaded: {self.api_key[:6]}...")
        self.base_url = "https://generativelanguage.googleapis.com/v1beta"
        self.model = "gemini-2.5-flash"
        
        logger.info(f"StudioService initialized with {self.model} (REST API v1beta)")

    def _extract_json(self, text: str) -> dict:
        """Safely extract JSON from response text"""
        start = text.find("{")
        end = text.rfind("}")
        if start == -1 or end == -1:
            raise ValueError("No JSON object found in response")
        return json.loads(text[start:end+1])

    async def process_command(self, prompt: str, context: Dict[str, Any], mode: str) -> Dict[str, Any]:
        """
        Processes a natural language command into a structured action.
        """
        if not self.api_key:
            return {
                "type": "NO_OP",
                "reason": "Gemini API key is missing. Please configure GEMINI_KEY in the backend environment."
            }

        try:
            # Prepare context-aware prompt
            context_str = json.dumps(context)
            full_prompt = f"{STUDIO_SYSTEM_PROMPT}\n\nMode: {mode}\nContext: {context_str}\n\nUser: {prompt}"
            
            # Call Gemini REST API directly
            url = f"{self.base_url}/models/{self.model}:generateContent"
            
            headers = {
                "Content-Type": "application/json",
            }
            
            payload = {
                "contents": [{
                    "role": "user",
                    "parts": [{
                        "text": full_prompt
                    }]
                }],
                "generationConfig": {
                    "temperature": 0.2,  # Low temperature for scientific accuracy
                    "topK": 40,
                    "topP": 0.95,
                    "maxOutputTokens": 1024,
                }
            }
            
            async with httpx.AsyncClient(timeout=30.0) as client:
                response = await client.post(
                    url,
                    params={"key": self.api_key},
                    headers=headers,
                    json=payload
                )
                
                if response.status_code != 200:
                    error_detail = response.text
                    logger.error(f"Gemini API error {response.status_code}: {error_detail}")
                    raise Exception(f"API returned {response.status_code}: {error_detail}")
                
                result = response.json()
                
                # Extract text from response
                if "candidates" in result and len(result["candidates"]) > 0:
                    candidate = result["candidates"][0]
                    if "content" in candidate and "parts" in candidate["content"]:
                        text = candidate["content"]["parts"][0]["text"]
                        
                        # Use safer JSON extraction
                        return self._extract_json(text)
                
                raise Exception("No valid response from Gemini")
                
        except Exception as e:
            # Log full error details for debugging
            logger.error(f"StudioService Error: {type(e).__name__}: {e}", exc_info=True)
            
            # Sanitize error for user
            error_msg = str(e)
            if "API_KEY_INVALID" in error_msg or "401" in error_msg or "403" in error_msg:
                user_msg = "AI service authentication failed. Please check API key configuration."
            elif "404" in error_msg or "not found" in error_msg:
                user_msg = "AI model endpoint error. Service may be temporarily unavailable."
            elif "429" in error_msg or "quota" in error_msg.lower():
                user_msg = "AI service rate limit reached. Please try again in a moment."
            else:
                user_msg = "AI service encountered an error. Please try again."
            
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
