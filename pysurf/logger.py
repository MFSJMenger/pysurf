from datetime import datetime
from functools import partial
#
from .utils.context_utils import DoOnException
from .utils.context_utils import ExitOnException
from .utils.context_utils import BaseContextDecorator


__all__ = ["get_logger", "Logger"]


def get_logger(filename, name, sublogger=None, mode="w"):
    """initialize a new logger"""
    fhandle = Logger.get_fhandle(filename, mode)
    return Logger(name, fhandle, sublogger)


class LogBlock(BaseContextDecorator):

    def __init__(self, logger, txt=None):
        self.txt = txt
        self.logger = logger

    def set_text(self, txt):
        self.txt = txt

    def __enter__(self):
        self.logger.info(f"\nEnter '{self.txt}' at: "
                         f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    def __exit__(self, exception_type, exception_value, traceback):
        self.logger.info(f"Leave '{self.txt}' at: "
                         f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")


class _LoggerBase(object):

    def __init__(self, name, fhandle):
        self.name = name
        self.fhandle = fhandle
        self._info_block = LogBlock(self)
        self._error_block = DoOnException(self.error)

    @staticmethod
    def get_fhandle(filename, mode="w"):
        """Generate a new file handle"""
        return _TextHandle(filename, mode)

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
        self.fhandle = handle

    def debug(self, txt):
        self.fhandle.write(f"Debug: {txt}\n")

    def header(self, name, dct=None):
        txt = '********************************************************************************\n'
        txt += '*{:^78}*\n'.format(name)
        txt += '*{:^78}*\n'.format(' ')
        txt += f"* {'Date':15}: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}".ljust(79) + '*\n'
        txt += '*{:^78}*\n'.format(' ')
        if dct is not None:
            for key in dct:
                inter = f"* {key:15}: {dct[key]}"
                if len(inter) < 80:
                    inter = inter.ljust(79) + '*\n'
                else:
                    inter_new = inter[0:79] + '*\n'
                    inter = inter[79:]
                    while len(inter) > 60:
                        inter_new += "* {:^15}  ".format(' ') + inter[0:60] + '*\n'
                        inter = inter[60:]
                    inter_new += "* {:^15}  ".format(' ') + inter.ljust(60) + '*\n'
                    inter = inter_new
                txt += inter
        txt += '*{:^78}*\n'.format(' ')
        txt += '********************************************************************************\n\n\n'
        self.fhandle.write(f"{txt}\n")


    def info(self, txt):
        self.fhandle.write(f"{txt}\n")

    def warning(self, txt):
        self.fhandle.write(f"Warning: {txt}\n")

    def info_block(self, txt):
        self._info_block.set_text(txt)
        return self._info_block

    def exit_on_exception(self, txt):
        self._error_block.set_args(txt)
        return self._error_block

    def error(self, txt):
        error = (f"in '{self.name}' at "
                 f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}:\n"
                 f"{txt}\n\n")
        #
        self.fhandle.write(f"Error Termination {error}")
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
            return {logger_name: Logger(f"{self.name.upper()}-{logger_name}",
                                        self.fhandle, sub_logger)
                    for logger_name, sub_logger in sublogger.items()}
        if isinstance(sublogger, tuple):
            return {logger: Logger(f"{self.name.upper()}-{logger}", self.fhandle)
                    for logger in sublogger}
        raise Exception("Sublogger can only be tuple, dict, list or None!")


class _TextHandle(object):

    def __init__(self, filename=None, mode="w"):
        self._setup(filename, mode=mode)

    def _setup(self, filename, mode="w"):
        if filename is None:
            self._f = None
            self.write = partial(print, end='')
        else:
            self._f = open(filename, mode=mode)
            self.write = self._write
    
    def _write(self, txt):
        self._f.write(txt)
        self._f.flush()

    def __del__(self):
        if self._f is not None:
            self._f.close()
