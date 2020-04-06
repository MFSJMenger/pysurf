import pytest

from pysurf.utils.context_utils import SetOnException, DoOnException, ExitOnException


@pytest.fixture
def dct():
    return {
            'x': 5,
            'y': 10,
            'z': [1, 2, 3],
            'k': 'hallo'
    }


def set_dct_value(dct, key, value):
    dct[key] = value


def test_do_on_exception_succeed(dct):

    with DoOnException(set_dct_value, dct, 'x', 100):
        set_dct_value(dct, 'x', 2)

    assert dct['x'] == 2


def test_do_on_exception_fail(dct):

    with DoOnException(set_dct_value, dct, 'x', 1111):
        raise Exception("")
        set_dct_value(dct, 'x', 2)

    assert dct['x'] == 1111


def test_set_on_exception_fail(dct):

    with SetOnException(dct) as f:
        f.set_value('x', 100)
        raise Exception("")
        f.set_value('y', [1, 2, 3])
        f.set_value('z', 'hallo')
        f.set_value('k', dict(a=10))
    x, y, z, k = f.result
    assert x == dct['x']
    assert y == dct['y']
    assert z == dct['z']
    assert k == dct['k']


def test_set_on_exception_succeed(dct):

    with SetOnException(dct) as f:
        f.set_value('x', 100)
        f.set_value('y', [1, 2, 3])
        f.set_value('z', 'hallo')
        f.set_value('k', dict(a=10))
    x, y, z, k = f.result
    assert x == 100
    assert y == [1, 2, 3]
    assert z == 'hallo'
    assert k == {'a': 10}


def test_exit_on_exception_fail():
    value = ""
    with pytest.raises(SystemExit) as pytest_error:
        with ExitOnException():
            raise Exception(value)
    assert pytest_error.type == SystemExit


def test_exit_on_exception_succeed():
    value = ""
    with ExitOnException():
        value = 25
    assert value == 25
