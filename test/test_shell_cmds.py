import pytest

from app.utils.shell_cmds import shell, stoperr


def test_shell():
    nonsense = "well hello there fine gentlemen, my name is batman"
    out = shell(f"echo '{nonsense}'", ret_output=True)
    assert out == b'well hello there fine gentlemen, my name is batman\n'


def test_stoperr():
    with pytest.raises(SystemError):
        stoperr("test")
