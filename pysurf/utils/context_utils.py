import sys
from functools import wraps
from collections import namedtuple


class BaseContextDecorator(object):

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

        

class DoOnException(BaseContextDecorator):
    """Performs fallback function, if exception is raised"""

    def __init__(self, fallback, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs
        self._fallback = fallback

    def set_args(self, *args, **kwargs):
        self._args = args
        self._kwargs = kwargs

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            self._fallback(*self._args, **self._kwargs)
            return True

class SetOnException(BaseContextDecorator):

    def __init__(self, dct, reset_all=True):
        """Context Manager to set defaults on exception,
           
           Args:
              dct, dict = dictionary containing all varibales and their corresponding
                          default values
        """
        self._tuple = namedtuple("RESULT", (key for key in dct))
        self._dct = {key: None for key in dct}
        self._defaults = dct
        self.result = None
        
    def set_value(self, name, value):
        if name in self._dct:
            self._dct[name] = value

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            self.result = self._tuple(*(value for value in self._defaults.values()))
            return True
        self.result = self._tuple(*(value if value is not None else self._defaults[key] for key, value in self._dct.items()))

    def __call__(self, func):
        raise NotImplementedError("SetOnException cannot be used as decorator")

class ExitOnException(BaseContextDecorator):

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            print(f"Error Termination: {exception_value}")
            sys.exit()
