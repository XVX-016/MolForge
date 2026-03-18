"""
Backend service package.

Keep this package import-light so route modules can import individual
services without eagerly loading optional ML / chemistry dependencies.
"""

__all__ = [
    "PredictionService",
    "MoleculeService",
    "GenerationService",
    "UserService",
]
