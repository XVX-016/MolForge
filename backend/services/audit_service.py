from datetime import datetime
from sqlmodel import Session
from typing import Dict, Any, Optional
from backend.models.db.audit import UserActionLog
import uuid

def log_user_action(
    db: Session,
    action_type: str,
    entity_type: str,
    entity_id: str,
    details: Dict[str, Any] = None,
    user_id: str = "anonymous",
    client_ip: str = None
):
    """
    Records an immutable audit log entry for regulatory compliance.
    """
    log_entry = UserActionLog(
        action_type=action_type,
        entity_type=entity_type,
        entity_id=str(entity_id),
        details=details or {},
        user_id=user_id,
        client_ip=client_ip
    )
    db.add(log_entry)
    # We commit immediately to ensure the log is persisted even if subsequent ops fail?
    # Or should we let the caller handle commit? 
    # For audit, we typically want it part of the transaction or flushed immediately.
    # Let's assume the caller manages the transaction scope for atomicity, 
    # OR we add it to the session and let the global commit handle it.
    # But for specialized audit logging, sometimes we want to guarantee logging independent of success if we are logging "attempts".
    # For "actions taken", it usually accompanies the state change.
    # Let's just add to DB session.
    # Note: caller must invoke db.commit()
    return log_entry
