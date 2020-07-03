from pysurf.workflow import engine


workflow = engine.create_workflow("sp_calc", """
crd = read_xyzfile_crd(crd_file)
atomids = read_xyzfile_atomids(crd_file)
spp = spp_calc("spp.inp", atomids, nstates, properties=properties)
res = sp_calc(spp, crd, properties=properties)
""")

wf = workflow.run()
#wf = workflow.run({"properties": ['energy']})

#print(wf['res']['energy'])
