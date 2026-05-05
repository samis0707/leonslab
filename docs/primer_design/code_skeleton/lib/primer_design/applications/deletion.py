"""Orchestrate the pEXG2 in-frame deletion pipeline.

Per skill_v2 §4 and project_plan §4.1.
"""
from __future__ import annotations

from .. import gene_finder, primer_designer, storage_adapter, vectors, verification
from ..exceptions import OffTargetDetected
from ..plasmid_builder import assemble
from ..types import DesignRequest, DesignResult


def run(request: DesignRequest) -> DesignResult:
    """End-to-end deletion design.

    Steps:
        1. Load gene record (bundled, no R2 fetch).
        2. Load vector + tail convention.
        3. Search primer set (exhaustive (N, C) × P1 × P2 × P3 × P4).
        4. Lazy-fetch genome from R2 for off-target scan.
        5. Off-target scan; raise on unintended product.
        6. Assemble final plasmid.
        7. Verify (hard checks).
        8. Build DesignResult and return.

    Phase-2 implementation: glue together the imported helpers. Each helper has
    its own pseudocode in the corresponding module + skill_v2 reference.
    """
    raise NotImplementedError(
        "Phase 2 step 7: orchestrate per skill_v2 §4 and tests/integration/"
        "test_deletion_LB001_lasB.py expected values."
    )
