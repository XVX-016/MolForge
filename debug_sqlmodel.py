import sys
import traceback
try:
    from backend.chemistry.models import Experiment, WorkflowNode, WorkflowNodeResult
    print("Imports successful")
    from sqlmodel import SQLModel, create_engine
    engine = create_engine("sqlite:///:memory:")
    SQLModel.metadata.create_all(engine)
    print("Table creation successful")
except Exception as e:
    print(f"Error type: {type(e).__name__}")
    print(f"Error message: {e}")
    traceback.print_exc()
    sys.exit(1)
