from pysurf.workflow import engine


workflow = engine.create_workflow("copy_execute", """
folders = get_subfolder(folder, subfolder)
copy_execute(folders, copy, exe)
""")

workflow.run()
