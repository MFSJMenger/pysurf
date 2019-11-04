from datetime import datetime
from functools import partial
import sys


__all__ = ["get_logger", "Logger"]


def get_logger(filename, name, handles):
    fhandle = _TextHandle(filename)
    return Logger(name, fhandle, handles)


class _LoggerBase(object):

    def __init__(self, name, fhandle):
        self.name = name
        self._enter = False
        self.fhandle = fhandle

    def set_fhandle(self, handle):
        """ """
        self.fhandle = _TextHandle(handle)

    @property
    def enter(self):
        self._enter = True
        self.enter_text(f"routine {self.name}")

    @property
    def leave(self):
        self._enter = False
        self.leave_text(f"routine {self.name}")

    def enter_text(self, txt):
        self.fhandle.write(f"\nEnter '{txt}' at: "
                           f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    def leave_text(self, txt):
        self.fhandle.write(f"\nLeave '{txt}' at: "
                           f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    def info(self, txt):
        if self._enter is True:
            txt = f"{txt}\n"
        else:
            txt = f"{self.name}: {txt}\n"
        self.fhandle.write(txt)

    def warning(self, txt):
        if self._enter is True:
            txt = f"Warning: {txt}\n"
        else:
            txt = f"{self.name} Warning: {txt}\n"
        self.fhandle.write(txt)

    def error(self, txt):
        error = (f"\n\n{self.name} Error Termination:\n"
                 f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                 f"{txt}\n\n")
        #
        self.fhandle.write(error)
        # raise Exception
        raise Exception(error)


class Logger(_LoggerBase):

    def __init__(self, name, fhandle, handles=None):
        _LoggerBase.__init__(self, name, fhandle)
        self.handles = self._get_handles(handles)

    def __getitem__(self, key):
        return self.handles[key]

    def add_handle(self, handles):
        """add new subhandles to logger"""
        if isinstance(handles, str):
            handles = [handles]
        #
        handels = self._get_handles(handles)
        # register new keys
        for key, value in handels.items():
            self.handles[ky] = value

    def _get_handles(self, handles):
        if handles is None:
            return
        if isinstance(handles, list):
            return {handle: Logger(f"{self.name.upper()}-{handle}", self.fhandle)
                    for handle in handles}
        if isinstance(handles, dict):
            return {key: Logger(f"{self.name.upper()}-{key}", self.fhandle, handle)
                    for key, handle in handles.items()}
        if isinstance(handles, tuple):
            return {handle: Logger(f"{self.name.upper()}-{handle}", self.fhandle)
                    for handle in handles}
        raise Exception("can only handle tuple, dict, list and None!")


class _TextHandle(object):

    def __init__(self, filename=None):
        self._setup(filename)

    def _setup(self, filename):
        if filename is None:
            self._f = None
            self.write = partial(print, end='')
        else:
            self._f = open(filename, "w")
            self.write = self._f.write

    def __del__(self):
        if self._f is not None:
            self._f.close()
