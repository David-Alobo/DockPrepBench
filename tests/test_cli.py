import pytest

from dockprepbench import __version__
from dockprepbench.cli import main


def test_main_returns_zero():
    assert main([]) == 0


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])

    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
