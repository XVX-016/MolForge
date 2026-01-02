import google.generativeai as genai
import os
import sys
from pathlib import Path

# Add backend directory to path
root = Path.cwd()
sys.path.append(str(root))
sys.path.append(str(root / "backend"))

from backend.config import settings

import asyncio

async def main():
    genai.configure(api_key=settings.GEMINI_KEY)
    model = genai.GenerativeModel(
        'models/gemini-1.5-flash',
        system_instruction="Output JSON."
    )

    print(f"Testing Gemini 1.5 Flash SYNC with JSON...")
    try:
        response = model.generate_content(
            "List 3 elements in JSON format with 'name' and 'symbol' fields.",
            generation_config={"response_mime_type": "application/json"}
        )
        print("SUCCESS SYNC JSON!")
        print(response.text)
    except Exception as e:
        print(f"FAILED SYNC JSON: {type(e).__name__}: {e}")

    print(f"Testing Gemini 1.5 Flash ASYNC...")
    try:
        response = await model.generate_content_async("Hi")
        print("SUCCESS ASYNC!")
    except Exception as e:
        print(f"FAILED ASYNC: {type(e).__name__}: {e}")

if __name__ == "__main__":
    asyncio.run(main())
