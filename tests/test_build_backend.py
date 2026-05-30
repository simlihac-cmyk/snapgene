import pytest

import plasmidlab_build_backend
from plasmidlab_build_backend import _ensure_supported_python


def test_build_backend_guard_accepts_python_312_plus() -> None:
    _ensure_supported_python((3, 12, 0))
    _ensure_supported_python((3, 13, 0))


def test_build_backend_guard_rejects_python_311() -> None:
    with pytest.raises(RuntimeError, match="requires Python 3.12 or newer"):
        _ensure_supported_python((3, 11, 9))


def test_build_backend_does_not_advertise_unsupported_editable_metadata_hook() -> None:
    assert not hasattr(plasmidlab_build_backend, "prepare_metadata_for_build_editable")
