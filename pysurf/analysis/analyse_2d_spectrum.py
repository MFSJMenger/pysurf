import numpy as np

from pysurf.workflow import engine
from pysurf.database import PySurfDB


@engine.register_action
def calc_2d_spec(files: "list", energy_start: "float", energy_end: "float", en_points: "int", timesteps: "int") -> "meshplotdata":
    result = np.zeros((en_points, timesteps))
    for file in files:
        db = PySurfDB.load_database(file, read_only=True)
        energy = np.copy(db['energy'])
        fosc = np.copy(db['fosc'])
        currstate = np.copy(db['currstate']).flatten()
        for idx, en in enumerate(energy):
            state = int(currstate[idx])
            en_diff = en[state] - en[0]
            en_pos = (en_diff-energy_start)/(energy_end-energy_start)*en_points
            en_pos = int(en_pos)
            if en_pos >=0 and en_pos < en_points:
                #units are arbitrary, therefor not full Einstain A coefficient ist used, but just f*v**2
                result[en_pos, idx] += fosc[idx][state]*en_diff**2
    X = np.arange(timesteps)*0.5
    Y = np.linspace(energy_start, energy_end, en_points)*27.2114
    return (X, Y, result)


workflow = engine.create_workflow("copy_execute", """
files = get_files(folder, subfolder, filename)
spec = calc_2d_spec(files, 0.075, 0.2, 100, 201)
plot = mesh_plot("plot_spec.inp", spec)
show_plot(plot)
""")

workflow.run()

