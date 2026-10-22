import numpy as np
import pytest

from dockprepbench.engines.base import Box, DockingResult


def test_docking_result_requires_one_score_per_pose():
    with pytest.raises(ValueError):
        DockingResult(poses=[], scores=[-5.0])


def test_box_coerces_numpy_scalars_to_float():
    box = Box(center=tuple(np.float32([1.5, 2, 3])), size=(np.float64(20), 20, 20))

    assert all(type(v) is float for v in box.center + box.size)
