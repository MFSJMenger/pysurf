from pysurf.colt import Colt
from pysurf.colt import FromCommandline
from pysurf.sh.run_trajectory import RunTrajectory




@FromCommandline("""
inputfile = propagation.inp :: file
""")
def command_run_trajectory(inputfile):
    RunTrajectory(inputfile)

if __name__=="__main__":
    command_run_trajectory()
