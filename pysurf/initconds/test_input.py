from pysurf.colt import AskQuestions

def ask_questions(inputfile):
        question_string = """
        # database file where the initial conditions are saved or from which the initial conditions
        # are taken if the file already exists.
        outputfile = initconds.db

        # Input source for the normal modes and/or frequencies, which are used to generate the 
        # initial conditions.
        # Possible options are:
        # - molden
        # - frequencies
        from = none :: str :: [molden, frequencies]

        # Describes which sampling algorithm is used to generate the initial conditions.
        # The default is wigner.
        sampling = wigner

        # Number of initial conditions that have to be created.
        # The default value is 100.
        number of initial conditions = 100 :: int

        # If initial conditions are generated from a molden file, this subquestion asks for the 
        # molden file.
        [from(molden)]
        moldenfile = none

        # If initial conditions are generated from frequencies, here a list of frequencies has to
        # be given for the modes. The list can be in the python list format or just comma separated
        # values.
        [from(frequencies)]
        frequencies = none
        """
    
        quests = AskQuestions.questions_from_string("INITIAL CONDITIONS", question_string, config=inputfile)
        print(quests['sampling'])
        

if __name__=="__main__":
    ask_questions('test.inp')
