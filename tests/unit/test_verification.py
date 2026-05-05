"""Tests for hard verification gates (D7.4)."""
from __future__ import annotations

import pytest

from primer_design.exceptions import VerificationFailure
# Detailed verification tests live in tests/integration/ where we have full
# DesignResult fixtures. Here we test the helpers in isolation when sensible.


class TestRecognitionCounting:
    """The recognition-count helper is private, but we exercise the whole
    verify_deletion path in tests/integration/test_deletion_LB001_lasB.py.
    Stub here for Phase-2 unit-level coverage of edge cases."""

    @pytest.mark.skip(reason="Phase 2: implement once verify() exposes a counted helper or factor it out")
    def test_palindromic_motif_not_double_counted(self):
        pass

    @pytest.mark.skip(reason="Phase 2")
    def test_zero_count_when_destroyed(self):
        pass

    @pytest.mark.skip(reason="Phase 2")
    def test_one_count_when_partial(self):
        pass
