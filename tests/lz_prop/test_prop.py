import numpy as np

from pysurf.dynamics import LandauZener
from pysurf.spp.model import HarmonicOscillator
from pysurf.sampling import DynCondition



class FakeDB():
    def __init__(self):
        self.crd = []
        self.veloc = []

    def add_step(self, time, data, v, iactive, ekin, epot, etot):
        self.crd.append(list(data['crd']))
        self.veloc.append(list(v))

    @property
    def model(self):
        return True

class FakeSPP:
    def __init__(self):
        self.ho = HarmonicOscillator()
        self.model = True

    def request(self, crd, properties, states):
        request = {'crd': crd}
        res = self.ho.get(request)
        return res

class LZ_Test(LandauZener):

    
    def __init__(self):
        self.spp = FakeSPP()
        self.restart = False
        self.masses = np.array([1.0])
        self.nstates = 1
        
        self.ekin_save = []
        self.epot_save = []
        self.etot_save = []
        
        self.init = DynCondition(crd=np.array([1.0]), veloc=np.array([0.0]), state=0)
        self.db = FakeDB()
        self._run(100, 0.1)

    def get_results(self):
        return np.array(self.db.crd).flatten(), np.array(self.db.veloc).flatten(), np.array(self.etot_save).flatten()

    def output_step(self, istep, time, iactive, ekin, epot, etot, diff):
        self.ekin_save.append(list([ekin]))
        self.epot_save.append(list([epot]))
        self.etot_save.append(list([etot]))

def test_lz_prop():
    crd, veloc, etot = LZ_Test().get_results()
    time = np.linspace(0, 10, 101)
    sin = -np.sin(time)
    cos = np.cos(time)
    
    assert(np.max(np.abs(crd - cos)) < 0.01)
    assert(np.max(np.abs(veloc - sin)) < 0.01)
    assert(np.max(np.abs(etot-0.5)))
        
