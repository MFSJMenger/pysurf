from pysurf.wigner import InitialConditions

conditions = InitialConditions.from_db("init_conds.db")

for i, con in enumerate(conditions):
    print(con)
    if i == 5:
        break

