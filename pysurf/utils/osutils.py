import os


def exists_and_isfile(filename):
    """if file does not exist, create it
       if file exisits check if isfile or  raise Exception
    """

    if os.path.exists(filename):
        if os.path.isfile(filename):
            return True
        else:
            raise Exception("Object '%s' exisits but is not a file!")

    return False
