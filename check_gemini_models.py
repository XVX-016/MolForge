
import os
import httpx
from pathlib import Path
from dotenv import load_dotenv

# Load env
env_path = Path('backend/.env')
load_dotenv(dotenv_path=env_path)

api_key = os.getenv('GEMINI_KEY') or os.getenv('GOOGLE_API_KEY')

if not api_key:
    print("ERROR: No API key found")
    exit(1)

def get_models(version):
    url = f"https://generativelanguage.googleapis.com/{version}/models?key={api_key}"
    try:
        response = httpx.get(url)
        if response.status_code == 200:
            return response.json().get('models', [])
        else:
            print(f"ERROR {version} ({response.status_code}): {response.text[:100]}")
            return []
    except Exception as e:
        print(f"EXCEPTION {version}: {e}")
        return []

print("=== Gemini API Model Check ===")
for ver in ['v1', 'v1beta']:
    print(f"\nAPI Version: {ver}")
    models = get_models(ver)
    for m in models:
        name = m.get('name', 'unknown')
        methods = m.get('supportedGenerationMethods', [])
        if 'generateContent' in methods:
            print(f"  [YES] {name}")
        else:
            print(f"  [NO ] {name} (Methods: {methods[:2]})")
