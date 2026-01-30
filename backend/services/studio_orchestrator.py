import asyncio
import logging
import uuid
from typing import Dict, Any
from sqlmodel import Session
from backend.db import engine
from backend.chemistry.models import WorkflowNode, WorkflowNodeResult, WorkflowNodeState

logger = logging.getLogger(__name__)

async def dispatch_node_execution(node_id: uuid.UUID, result_id: uuid.UUID):
    """
    Background orchestrator that simulates scientific execution.
    """
    logger.info(f"Starting execution of Node {node_id} (Result {result_id})")
    
    with Session(engine) as db:
        result = db.get(WorkflowNodeResult, result_id)
        node = db.get(WorkflowNode, node_id)
        if not result or not node:
            logger.error("Node or Result not found in background task")
            return

        # 1. Update to RUNNING
        result.status = WorkflowNodeState.RUNNING
        db.add(result)
        db.commit()
        
        try:
            # 2. Simulate scientific calculation
            # Each node type has its own execution logic
            if node.node_type == "docking":
                await simulate_docking(node, result)
            elif node.node_type == "md":
                await simulate_md(node, result)
            else:
                await simulate_generic(node, result)
            
            # 3. Mark as COMPLETED
            result.status = WorkflowNodeState.COMPLETED
            logger.info(f"Node {node_id} completed successfully")
            
        except Exception as e:
            logger.error(f"Node {node_id} failed: {e}")
            result.status = WorkflowNodeState.FAILED
            result.stdout_logs = str(e)
        
        db.add(result)
        db.commit()

async def simulate_docking(node: WorkflowNode, result: WorkflowNodeResult):
    """Simulate a docking run"""
    await asyncio.sleep(5) # Simulate compute time
    result.output_data = {
        "binding_affinity": -8.4,
        "poses_url": "s3://molforge/poses/pose_1.sdf",
        "interaction_residues": ["HIS12", "TRP45"],
        "engine": "AutoDock Vina V1.2"
    }
    result.artifacts_checksum = "sha256:abc123docking"

async def simulate_md(node: WorkflowNode, result: WorkflowNodeResult):
    """Simulate a Molecular Dynamics run"""
    await asyncio.sleep(10) # MD takes longer
    result.output_data = {
        "rmsd_url": "s3://molforge/plots/rmsd.csv",
        "stability_score": 0.88,
        "trajectory_url": "s3://molforge/trajectories/traj.dcd"
    }
    result.artifacts_checksum = "sha256:xyz789md"

async def simulate_generic(node: WorkflowNode, result: WorkflowNodeResult):
    """Fallback for unknown node types"""
    await asyncio.sleep(2)
    result.output_data = {"status": "success", "message": "Generic completion"}
