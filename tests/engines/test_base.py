import pytest

from dockprepbench.engines.base import DockingResult


def test_docking_result_requires_one_score_per_pose():
    with pytest.raises(ValueError):
        DockingResult(poses=[], scores=[-5.0])
