FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    cmake \
    git \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"

# Accept Conda Terms of Service for main and r channels
RUN conda config --set always_yes yes --set changeps1 no && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install RDKit + Python 3.11 from conda-forge
RUN conda install -y -c conda-forge rdkit=2023.09 python=3.11 && \
    conda clean -afy

# Copy requirements (if any)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy app code
COPY . .

EXPOSE $PORT
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
