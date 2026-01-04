
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

def check_version(version):
    print(f"\n--- Checking API version: {version} ---")
    url = f"https://generativelanguage.googleapis.com/{version}/models?key={api_key}"
    try:
        response = httpx.get(url)
        if response.status_code == 200:
            models = response.json().get('models', [])
            for m in models:
                name = m.get('name')
                methods = m.get('supportedGenerationMethods', [])
                if 'generateContent' in methods:
                    print(f"  [OK]  {name}")
                else:
                    print(f"  [--]  {name} (Methods: {methods})")
        else:
            print(f"  ERROR {response.status_code}: {response.text}")
    except Exception as e:
        print(f"  EXCEPTION: {e}")

check_version('v1')
check_version('v1beta')
