"""
Studio Service
Molecular design orchestrator powered by Gemini 1.5 Pro
"""
import os
import logging
import json
import re
from typing import Dict, Any, Optional
from rdkit import Chem
from backend.services.mentor_service import GEMINI_AVAILABLE

if GEMINI_AVAILABLE:
    import google.generativeai as genai

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
        if not GEMINI_AVAILABLE:
            raise ValueError("google-generativeai package not installed")
        
        api_key = os.getenv("GEMINI_KEY")
        if not api_key:
            raise ValueError("GEMINI_KEY not found in environment variables")
        
        genai.configure(api_key=api_key)
        self.model = genai.GenerativeModel(
            model_name="gemini-1.5-pro",
            system_instruction=STUDIO_SYSTEM_PROMPT
        )
        logger.info("StudioService initialized")

    async def process_command(self, prompt: str, context: Dict[str, Any], mode: str) -> Dict[str, Any]:
        """
        Processes a natural language command into a structured action.
        """
        try:
            # Prepare context-aware prompt
            context_str = json.dumps(context)
            full_prompt = f"Mode: {mode}\nContext: {context_str}\n\nUser: {prompt}"
            
            response = self.model.generate_content(full_prompt)
            # Try to extract JSON from markdown if Gemini wraps it
            text = response.text
            json_match = re.search(r"\{.*\}", text, re.DOTALL)
            if json_match:
                return json.loads(json_match.group())
            
            return json.loads(text)
        except Exception as e:
            logger.error(f"StudioService Error: {e}", exc_info=True)
            return {
                "type": "NO_OP",
                "reason": f"AI Orchestrator failed: {str(e)}"
            }

_studio_service: Optional[StudioService] = None

def get_studio_service() -> StudioService:
    global _studio_service
    if _studio_service is None:
        _studio_service = StudioService()
    return _studio_service
