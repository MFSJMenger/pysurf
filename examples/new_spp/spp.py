from pysurf import SurfacePointProvider


spp = SurfacePointProvider('spp.inp', ['energy', 'gradient'], 2, 3, ['h', 'h']) 



print(spp.request([[1, 0, 0], [0, 0, 0]], ['gradient'], [2]))
