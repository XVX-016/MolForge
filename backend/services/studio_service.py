
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
You are the MolForge AI Commander. Your role is strictly to PARSE intent from user requests into structured commands for the RDKit Chemistry Kernel.
You MUST ONLY output a JSON object. No markdown, no prose, no explanations outside the JSON.

COMMAND CONTRACT:
1. CREATE_MOLECULE: { "smiles": "SMILES_STRING" }
2. ADD_ATOM: { "element": "C", "position": [x, y, z] }
3. REPLACE_ATOM: { "atomId": "id", "newElement": "N" }
4. REMOVE_ATOM: { "atomId": "id" }
5. ADD_BOND: { "from": "id1", "to": "id2", "order": 1 }
6. REMOVE_BOND: { "bondId": "id" }
7. OPTIMIZE_GEOMETRY: {}
8. ANALYZE_PROPERTIES: {}
9. NO_OP: { "reason": "Scientific justification or error explanation" }

STRATEGIC RULES:
- You are an ORCHESTRATOR. You do not compute properties; RDKit does.
- If the user asks for "LogP" or "Stability", use ANALYZE_PROPERTIES.
- If the request is ambiguous, use NO_OP with a request for clarification in "reason".
- Never invent atomic coordinates unless performing a CREATE_MOLECULE from scratch.
- For structural edits, refer to the provided Molecule Context for atom/bond IDs.

OUTPUT FORMAT:
{
  "type": "ACTION_TYPE",
  "payload": { ... },
  "reason": "Technical justification for this design choice."
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
            logger.info("StudioService initialized with Gemini v1beta REST API")

    def _extract_json(self, response_data: dict) -> dict:
        """Robustly extract JSON from Gemini REST response"""
        try:
            parts = response_data.get('candidates', [{}])[0].get('content', {}).get('parts', [])
            if not parts:
                raise ValueError("No content returned from AI.")
            
            text = parts[0].get('text', '').strip()
            
            # Remove markdown code blocks if present
            if text.startswith("```json"):
                text = text[7:]
            if text.startswith("```"):
                text = text[3:]
            if text.endswith("```"):
                text = text[:-3]
            
            text = text.strip()
            
            # Find JSON block boundaries just in case
            start = text.find("{")
            end = text.rfind("}")
            if start != -1 and end != -1:
                return json.loads(text[start:end+1])
            
            return json.loads(text)
        except (KeyError, IndexError, json.JSONDecodeError) as e:
            logger.error(f"Failed to parse Gemini response: {e}. Raw text: {text[:200] if 'text' in locals() else 'N/A'}")
            return {
                "type": "NO_OP",
                "reason": f"AI Parsing Failure: Expected JSON command, received malformed output. Chemistry kernel standing by."
            }

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
            # Force use 2.0-flash-exp on v1beta which we verified works for this key
            model_name = "gemini-2.0-flash-exp"
            url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={self.api_key}"
            
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

            print(f"DEBUG: Calling Gemini API URL: {url}")
            print(f"DEBUG: Payload: {json.dumps(payload, indent=2)}")

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
