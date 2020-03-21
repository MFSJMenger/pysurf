import os
import numpy as np

from pysurf.database.database import Database
from pysurf.colt import FromCommandline

@FromCommandline("""
outfile = population.dat :: file
""")
def get_population_command(outfile):
    get_population(outfile)

def get_population(outfile):
    lsdir = os.listdir('./')
    population = [] 
    counter = []
    for obj in lsdir:
        if obj.startswith('traj.') and os.path.isfile(obj + '/prop.db'):
            try:
                db = Database.load_db(obj + '/prop.db')
                nstates = len(db['energy'][0])
                for step, st in enumerate(db['curr_state']):
                    pop = np.zeros(nstates, dtype=int)
                    pop[int(st)] = 1
                    if step >= len(population):
                        population += [pop]
                        counter += [1]
                    else:
                        population[step] = population[step] + pop
                        counter[step] += 1
            except:
                pass
    
    for i in range(len(population)):
        population[i] = population[i]/counter[i]
    
    with open(outfile, 'w') as output:
        for step, pop in enumerate(population):
            string = ''
            string += '{}   '.format(step)
            for p in pop:
                string += '{0:12.8f}   '.format(p)
            string += '\n'
            output.write(string)
       

if __name__=="__main__":
    get_population_command()
