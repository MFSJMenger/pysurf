import os


def exists_and_isfile(filename):
    """if file does not exist return False
       if file exisits check if isfile or raise Exception
    """
    if filename is None:
        return False
    if os.path.exists(filename):
        if os.path.isfile(filename):
            return True
        else:
            raise Exception("Object '%s' exisits but is not a file!")
    return False
