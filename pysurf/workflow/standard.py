from . import engine


@engine.register_action
def print(arg: 'anything'):
    myprint = globals()['__builtins__']['print']
    myprint(arg)
