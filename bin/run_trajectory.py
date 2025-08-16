from colt import from_commandline
from pysurf.fssh import State, VelocityVerletPropagator




@from_commandline("""
inputfile = prop.inp :: file
""")
def command_run_trajectory(inputfile="prop.inp"):
    elec_state = State.from_questions(config=inputfile)
    DY = VelocityVerletPropagator(elec_state)
    try:
        result_2 = DY.run()
    except SystemExit as err:
        print("An error:", err)


if __name__=="__main__":
    command_run_trajectory()
