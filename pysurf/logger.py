import sys
from datetime import datetime
from functools import partial
from .utils.context_utils import ExitOnException


__all__ = ["get_logger", "Logger"]


def get_logger(filename, name, sublogger):
    """initialize a new logger"""
    fhandle = Logger.get_fhandle(filename)
    return Logger(name, fhandle, sublogger)


class _LoggerBase(object):

    def __init__(self, name, fhandle):
        self.name = name
        self._enter = False
        self.fhandle = fhandle

    @staticmethod
    def get_fhandle(filename):
        """Generate a new file handle"""
        return _TextHandle(filename)

    def set_fhandle(self, handle):
        """set fhandle to new handle"""
        if isinstance(handle, _TextHandle):
            pass
        elif isinstance(handle, str) or isinstance(handle, None): 
            handle = self.get_fhandle(handle)
        else:
            raise Exception("new file handle can only be: \n"
                            "None     -> Console logger\n"
                            "filename -> File logger \n"
                            "logger   -> logger \n")

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

    def debug(self, txt):            
        self.fhandle.write(f"Debug: {txt}\n")

    def info(self, txt):
        self.fhandle.write(f"{txt}\n")

    def warning(self, txt):
        self.fhandle.write(f"Warning: {txt}\n")

    def error(self, txt):
        error = (f"\n\n{self.name} Error Termination:\n"
                 f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                 f"{txt}\n\n")
        #
        self.fhandle.write(error)
        # raise Exception
        with ExitOnException():
            raise Exception(error)


class Logger(_LoggerBase):

    def __init__(self, name, fhandle=None, handles=None):
        if fhandle is None:
            # create an terminal logger
            fhandle = self.get_fhandle(None)
        _LoggerBase.__init__(self, name, fhandle)
        self.sublogger = self._generate_sublogger(handles)

    def __getitem__(self, key):
        return self.sublogger[key]

    def add_sublogger(self, sublogger):
        """add new sublogger to logger"""
        if sublogger is None:
            return
        if isinstance(sublogger, str):
            sublogger = [sublogger]
        #
        sublogger = self._generate_sublogger(sublogger)
        # register new keys
        for key, value in sublogger.items():
            self.sublogger[key] = value

    def _generate_sublogger(self, sublogger):
        """create subloggers"""
        if sublogger is None:
            return
        if isinstance(sublogger, list):
            return {logger: Logger(f"{self.name.upper()}-{logger}", self.fhandle)
                    for logger in sublogger}
        if isinstance(sublogger, dict):
            return {logger_name: Logger(f"{self.name.upper()}-{logger_name}", self.fhandle, sub_logger)
                    for logger_name, sub_logger in sublogger.items()}
        if isinstance(sublogger, tuple):
            return {logger: Logger(f"{self.name.upper()}-{logger}", self.fhandle)
                    for logger in sublogger}
        raise Exception("Sublogger can only be tuple, dict, list or None!")


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
