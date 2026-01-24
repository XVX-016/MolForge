from typing import Dict, Any, Optional
import logging
from rdkit import rdBase

logger = logging.getLogger(__name__)

# Cache Versioning
RDKIT_VERSION = rdBase.rdkitVersion
PROPS_VERSION = "v1"

class DescriptorCache:
    def __init__(self, max_size: int = 5000):
        self._cache: Dict[str, Dict[str, Any]] = {}
        self.max_size = max_size

    def _make_key(self, inchikey: str) -> str:
        return f"{inchikey}:{RDKIT_VERSION}:{PROPS_VERSION}"

    def get(self, inchikey: str) -> Optional[Dict[str, Any]]:
        key = self._make_key(inchikey)
        return self._cache.get(key)

    def set(self, inchikey: str, properties: Dict[str, Any]):
        if len(self._cache) >= self.max_size:
            # Simple FIFO eviction
            first_key = next(iter(self._cache))
            del self._cache[first_key]
            
        key = self._make_key(inchikey)
        self._cache[key] = properties
        logger.debug(f"Cached properties for identifier: {key}")

# Singleton instance
descriptor_cache = DescriptorCache()
