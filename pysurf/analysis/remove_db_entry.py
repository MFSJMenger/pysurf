from pysurf.database.database import Database
from colt import from_commandline


@from_commandline("""
db = :: file
rm_entry = :: int
""")
def remove_entry_command(db, rm_entry)
    remove_entry(db, rm_entry)

def remove_entry(dbfile, rm_entry):
    rm_entry = int(rm_entry)
    print('database: ', dbfile)
    print('remove entry: ', rm_entry)

    db = Database.load_db(dbfile)
    dbnew = Database.empty_like('db_new.dat', db)
    
    keys = db.get_keys()
    for key in keys:
        lendb = len(db[key])

    for i in range(lendb):
        if i == rm_entry:
            print(f'Do not copy entry: {i}')
            continue

        for key in keys:
            dbnew.append(key, db[key][i])
        dbnew.increase

if __name__=='__main__':
    remove_entry_command()
