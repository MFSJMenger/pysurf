from . import engine

@engine.register_action(["style"])
def print_wf(anything):
    print(anything)
