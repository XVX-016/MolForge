# Optimized Dockerfile for Railway - Fast builds
# Build context MUST be set to /backend in Railway Settings → Deployment → Directory
# This Dockerfile assumes Railway build context is already /backend

FROM python:3.11-slim

WORKDIR /app

# Install minimal build dependencies
# Note: rdkit-pypi comes with pre-compiled binaries, but some operations may need gcc/g++
# If you get RDKit import errors, uncomment Boost libraries below
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Optional: Uncomment if RDKit operations fail without Boost
# Most RDKit operations work fine without system Boost (rdkit-pypi includes it)
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     libboost-system1.74.0 \
#     libboost-thread1.74.0 \
#     && rm -rf /var/lib/apt/lists/* \
#     && apt-get clean

# Copy requirements first for better Docker layer caching
COPY requirements.txt .

# Upgrade pip and install Python dependencies
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy backend code (context is already /backend, so no backend/ prefix needed)
COPY . .

# Set Python path for imports
ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1

# Railway provides PORT environment variable
ENV PORT=8000

# Expose port
EXPOSE ${PORT}

# Health check (using curl instead of requests to avoid extra dependency)
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:${PORT:-8000}/health')" || exit 1

# Run the application
CMD ["sh", "-c", "uvicorn app:app --host 0.0.0.0 --port ${PORT:-8000}"]
