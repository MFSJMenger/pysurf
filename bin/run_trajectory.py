from pysurf.colt import Colt
from pysurf.colt import from_commandline
from pysurf.dynamics.run_trajectory import RunTrajectory




@from_commandline("""
inputfile = prop.inp :: file
""")
def command_run_trajectory(inputfile):
    RunTrajectory.from_inputfile(inputfile)

if __name__=="__main__":
    command_run_trajectory()
