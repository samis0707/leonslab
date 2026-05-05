"""Application dispatch."""
from __future__ import annotations

from . import deletion, expression, tagging
from ..types import DesignRequest, DesignResult


def dispatch(request: DesignRequest) -> DesignResult:
    """Route to the right application orchestrator based on ``request.application``."""
    app = request.application
    if app == "deletion":
        return deletion.run(request)
    if app == "tagging":
        return tagging.run(request)
    if app == "expression":
        return expression.run(request)
    raise ValueError(f"Unknown application: {app}")
