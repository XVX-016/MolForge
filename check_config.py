import os
import sys
from pathlib import Path

# Add backend directory to path
root = Path.cwd()
sys.path.append(str(root))
sys.path.append(str(root / "backend"))

from backend.config import settings

print(f"DEBUG SCRIPT: CWD is {os.getcwd()}")
key = settings.GEMINI_KEY
print(f"DEBUG SCRIPT: GEMINI_KEY exists? {bool(key)}")
if key:
    print(f"DEBUG SCRIPT: GEMINI_KEY length: {len(key)}")
    print(f"DEBUG SCRIPT: GEMINI_KEY repr: {repr(key)}")
    print(f"DEBUG SCRIPT: GEMINI_KEY starts with AIza? {key.startswith('AIza')}")
    print(f"DEBUG SCRIPT: GEMINI_KEY contains space/newline? {any(c.isspace() for c in key)}")

# Try a very simple generate call to confirm
import google.generativeai as genai
genai.configure(api_key=key)
try:
    model = genai.GenerativeModel('gemini-1.5-flash')
    resp = model.generate_content('Hi')
    print("DEBUG SCRIPT: Generation SUCCESS")
except Exception as e:
    print(f"DEBUG SCRIPT: Generation FAILED: {type(e).__name__}: {e}")
