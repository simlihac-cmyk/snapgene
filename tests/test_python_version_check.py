from scripts.check_python import is_supported, version_label


def test_python_version_check_accepts_supported_versions() -> None:
    assert is_supported((3, 12, 0))
    assert is_supported((3, 13, 0))


def test_python_version_check_rejects_unsupported_versions() -> None:
    assert not is_supported((3, 11, 9))
    assert not is_supported((2, 7, 18))


def test_python_version_label_uses_major_minor_micro() -> None:
    assert version_label((3, 12, 4, "final", 0)) == "3.12.4"
