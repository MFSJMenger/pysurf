from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.dynamics.run_trajectory import RunTrajectory




@FromCommandline("""
inputfile = prop.inp :: file
""")
def command_run_trajectory(inputfile):
    RunTrajectory.from_inputfile(inputfile)

if __name__=="__main__":
    command_run_trajectory()
