import numpy.random as random
from pysurf.wigner import WignerSampling
from pysurf.wigner import InitialConditions


# set your random number seed 
random.seed(1000)
# setup wigner sampling from molden file
sampling = WignerSampling.from_molden('molden.in')
# create an intial conditions file
NConds = 1000
conditions = sampling.create_initial_conditions('init_conds.db', NConds)
# 
mol = conditions.molecule

equilibrium = conditions.equilibrium
print(equilibrium)
# iterate over all conditions
for con in conditions:
    print(con)
    break
# get a specific intial condition
cond = conditions.get_condition(10)
# Read existing conditions
condition = InitialConditions.from_db("init_conds.db", 10, 24)
