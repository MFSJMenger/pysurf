import sys
from functools import wraps


class _BaseContextDecorator(object):

    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __call__(self, func):

        @wraps(func)
        def _wrapper(*args, **kwargs):
            with self:
                return func(*args, **kwargs)

        return _wrapper


class DoOnException(_BaseContextDecorator):
    """Performs fallback function, if exception is raised"""

    def __init__(self, fallback, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self._fallback = func

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            self._fallback(*self._args, **self._kwargs)
            return True


class ExitOnException(_BaseContextDecorator):

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            print(f"Error Termination: {exception_value}")
            sys.exit()
