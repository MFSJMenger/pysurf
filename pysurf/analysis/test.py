from colt import Colt


class Test(Colt):
    _questions = """
        testfloat = :: python(dict)
    """

    @classmethod
    def from_inputfile(cls, inputfile):
        quests = cls.generate_input(inputfile, config=None)
        quests = cls.generate_questions(config=inputfile)
        config = quests.check_only(inputfile)
        print(config['testfloat'][1])

if __name__=="__main__":
        Test.from_inputfile('test.inp')
