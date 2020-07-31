from shutil import copy2 as copy
from subprocess import run, CalledProcessError

from pysurf.utils import SubfolderHandle, FileHandle

from . import engine

@engine.register_action(need_self=True, iterator_id=0, progress_bar=True)
def copy_execute(self, folder: "list", filenames: "list", exe: "str"):
    for item in filenames:
        copy(item, folder)
    try:
        run(exe, cwd=folder, check=True, shell=True)
    except KeyboardInterrupt or CalledProcessError:
        self.error("Keyboard Interrupt")

@engine.register_action
def get_subfolder(folder: "str", subfolder: "str") -> "list":
    return SubfolderHandle(folder, subfolder)

@engine.register_action
def get_files(folder: "str", subfolder: "str", filename: "str") -> "list":
    return FileHandle(folder, subfolder, filename)


@engine.register_action
def test_subfolder(folder: "str", subfolder: "str"):
    pass
