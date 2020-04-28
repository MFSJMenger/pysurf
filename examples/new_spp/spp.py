from pysurf import SurfacePointProvider
from pysurf.utils import exists_and_isfile
from pysurf.colt import AskQuestions
b = AskQuestions(SurfacePointProvider.questions)
print(b.literals.data)
print(b.literals._literals)
SurfacePointProvider.generate_input('spp.inp', config=None)
spp = SurfacePointProvider('spp.inp', ['energy', 'gradient'], 2, 3, ['h', 'h']) 



print(spp.request([[1, 0, 0], [0, 0, 0]], ['gradient']))
