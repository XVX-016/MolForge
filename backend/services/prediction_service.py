# backend/services/prediction_service.py
"""
Prediction Service - High-level service for molecular property predictions
"""

from typing import Dict, List, Optional, Any
import logging
from backend.ml.prediction_engine import PredictionEngine

logger = logging.getLogger(__name__)


class PredictionService:
    """
    Service layer for molecular property prediction.
    Wraps the PredictionEngine with convenient methods for API usage.
    """
    
    _engine: Optional[PredictionEngine] = None
    
    @classmethod
    def _get_engine(cls) -> PredictionEngine:
        """Get singleton prediction engine instance."""
        if cls._engine is None:
            cls._engine = PredictionEngine()
        return cls._engine
    
    @classmethod
    def predict_properties(
        cls,
        smiles: str,
        properties: Optional[List[str]] = None,
        model_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Predict molecular properties from SMILES string.
        
        Args:
            smiles: SMILES string representation of the molecule
            properties: List of properties to predict (optional, defaults to all)
            model_id: Specific model to use (optional, uses default)
        
        Returns:
            Dictionary containing predictions and metadata
        """
        try:
            engine = cls._get_engine()
            
            # Prepare input data
            input_data = {"smiles": smiles}
            
            # Make prediction
            result = engine.predict(
                input_data=input_data,
                model_id=model_id,
                properties=properties,
            )
            
            # Convert to dictionary format expected by tasks
            return result.to_dict()
            
        except Exception as e:
            logger.error(f"Prediction failed for SMILES {smiles}: {e}")
            
            # Return fallback response with error indication
            fallback_props = properties or ["logP", "solubility", "toxicity", "molecularWeight"]
            predictions = {prop: 0.0 for prop in fallback_props}
            
            return {
                "predictions": predictions,
                "model_id": None,
                "warnings": [f"Prediction failed: {str(e)}"],
                "error": str(e)
            }
    
    @classmethod
    def predict_batch(
        cls,
        smiles_list: List[str],
        properties: Optional[List[str]] = None,
        model_id: Optional[str] = None,
        batch_size: int = 32,
    ) -> List[Dict[str, Any]]:
        """
        Batch predict properties for multiple SMILES.
        
        Args:
            smiles_list: List of SMILES strings
            properties: List of properties to predict
            model_id: Specific model to use
            batch_size: Batch size for processing
        
        Returns:
            List of prediction dictionaries
        """
        try:
            engine = cls._get_engine()
            
            # Prepare input data
            inputs = [{"smiles": smiles} for smiles in smiles_list]
            
            # Make batch predictions
            results = engine.predict_batch(
                inputs=inputs,
                model_id=model_id,
                batch_size=batch_size,
            )
            
            # Convert to dictionary format
            return [result.to_dict() for result in results]
            
        except Exception as e:
            logger.error(f"Batch prediction failed: {e}")
            
            # Return fallback responses
            fallback_props = properties or ["logP", "solubility", "toxicity", "molecularWeight"]
            fallback_predictions = {prop: 0.0 for prop in fallback_props}
            
            return [{
                "predictions": fallback_predictions,
                "model_id": None,
                "warnings": [f"Prediction failed: {str(e)}"],
                "error": str(e)
            } for _ in smiles_list]
    
    @classmethod
    def get_available_models(cls) -> List[Dict[str, Any]]:
        """
        Get list of available prediction models.
        
        Returns:
            List of model information dictionaries
        """
        try:
            engine = cls._get_engine()
            return engine.registry.list_models()
        except Exception as e:
            logger.error(f"Failed to get available models: {e}")
            return []
    
    @classmethod
    def get_available_properties(cls) -> List[str]:
        """
        Get list of available molecular properties that can be predicted.
        
        Returns:
            List of property names
        """
        # Default properties supported by the prediction engine
        return [
            "logP",           # Lipophilicity
            "solubility",     # Water solubility 
            "toxicity",       # Toxicity score
            "molecularWeight", # Molecular weight
        ]
    
    @classmethod
    def clear_cache(cls):
        """Clear prediction engine cache."""
        if cls._engine:
            cls._engine.clear_cache()
    
    @classmethod
    def reload_models(cls):
        """Reload prediction models (useful for model updates)."""
        if cls._engine:
            # Clear loaded models to force reload
            cls._engine.loaded_models.clear()
        cls._engine = None  # Force engine recreation
