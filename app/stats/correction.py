"""Multiple-testing correction over a vector of p-values."""

from __future__ import annotations

import numpy as np

from app.models import CorrectionMethod


def adjust_pvalues(p_values: np.ndarray, method: CorrectionMethod) -> np.ndarray:
    """Return adjusted p-values aligned with the input order."""
    raise NotImplementedError
