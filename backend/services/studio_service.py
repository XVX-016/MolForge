
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
You are the MolForge Studio Planner.

Your role is strictly limited to selecting and sequencing deterministic chemistry rules.
You do NOT generate molecules, SMILES, reactions, or chemical facts.

Rules:
- You may ONLY select from the provided rule set.
- You may ONLY reference properties provided in the analysis context.
- You MUST NOT invent chemical structures or transformations.
- You MUST NOT output free text outside the JSON schema.
- You MUST NOT suggest actions that violate HIGH severity alerts.

You act as a constrained planner that proposes optimization steps.
All chemistry execution is performed by the deterministic RDKit engine.
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

    async def process_command(self, prompt: str, context: Dict[str, Any], mode: str, analysis_context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Processes a natural language command using Gemini v1 REST API.
        Injects deterministic analysis to prevent hallucinations.
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
            analysis_str = json.dumps(analysis_context or {"issues": [], "suggestions": []})
            system_instruction = STUDIO_SYSTEM_PROMPT
            
            payload = {
                "contents": [
                    {
                        "role": "user",
                        "parts": [
                            {"text": f"SYSTEM_INSTRUCTION:\n{system_instruction}\n\nMODE: {mode}\nMOLECULE_CONTEXT: {context_str}\nANALYSIS_CONTEXT: {analysis_str}\n\nUSER_PROMPT: {prompt}"}
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
                action = self._extract_json(data)
                
                # Hard Validation: Enforce Allowed Commands & Reject Chemistry
                raw_text = json.dumps(action)
                if any(x in raw_text.upper() for x in ["C1=", "C=", "N1="]) or "SMILES" in raw_text.upper():
                    logger.error("AI attempted to inject SMILES/Chemistry truth. REJECTED.")
                    return {
                        "type": "NO_OP",
                        "reason": "AI attempted to invent chemistry truth (SMILES). Orchestration layer blocked this violation."
                    }

                allowed_types = ["select_rule", "propose_graph", "explain", "NO_OP"]
                if action.get("type") not in allowed_types:
                    logger.warning(f"AI attempted invalid command type: {action.get('type')}")
                    return {
                        "type": "NO_OP",
                        "reason": f"AI attempted an unauthorized command: {action.get('type')}. Rules engine requires strict JSON commands."
                    }
                
                # Propose Graph Validation
                if action.get("type") == "propose_graph":
                    steps = action.get("steps", [])
                    if len(steps) > 3:
                        return {
                            "type": "NO_OP",
                            "reason": "AI proposal graph exceeded complexity limits (max 3 steps). Refine your intent."
                        }

                # Rule Identification
                if action.get("type") == "select_rule":
                    rule_id = action.get("rule_id")
                    valid_ids = [s.get("id") for s in (analysis_context or {}).get("suggestions", [])]
                    if rule_id not in valid_ids:
                        logger.warning(f"AI hallucinated rule_id: {rule_id}")
                        return {
                            "type": "NO_OP",
                            "reason": f"AI rule '{rule_id}' rejected. Deterministic kernel only recognizes registered rules."
                        }
                
                return action
                
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
