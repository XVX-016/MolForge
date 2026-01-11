
import os
import httpx
import json
from pathlib import Path
from dotenv import load_dotenv

# Load env from backend/.env
env_path = Path('backend/.env')
load_dotenv(dotenv_path=env_path, override=True)

api_key = os.getenv('GEMINI_KEY') or os.getenv('GOOGLE_API_KEY')

if not api_key:
    print("ERROR: No API key found (Checked GEMINI_KEY and GOOGLE_API_KEY)")
    exit(1)

# Clean key
api_key = "".join(api_key.split())

def get_models(version):
    url = f"https://generativelanguage.googleapis.com/{version}/models?key={api_key}"
    try:
        response = httpx.get(url, timeout=10.0)
        if response.status_code == 200:
            return response.json().get('models', [])
        else:
            print(f"ERROR {version} ({response.status_code}): {response.text}")
            return []
    except Exception as e:
        print(f"EXCEPTION {version}: {e}")
        return []

print(f"=== Gemini API Model Check (Key: {api_key[:5]}...{api_key[-5:]}) ===")
for ver in ['v1', 'v1beta']:
    print(f"\nAPI Version: {ver}")
    models = get_models(ver)
    if not models:
        print("  No models found or error occurred.")
        continue
    
    for m in models:
        name = m.get('name', 'unknown')
        short_name = name.split('/')[-1]
        methods = m.get('supportedGenerationMethods', [])
        
        status = "[YES]" if 'generateContent' in methods else "[ NO]"
        print(f"  {status} {name} (Short: {short_name})")
        if 'generateContent' not in methods:
            print(f"       Methods: {methods}")

print("\n=== Recommendation ===")
print("If you see 'models/gemini-1.5-flash' listed with [YES], use it.")
print("The URL should be: https://generativelanguage.googleapis.com/{VERSION}/{NAME}:generateContent?key=...")
