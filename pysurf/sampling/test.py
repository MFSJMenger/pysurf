
class testsuper():
    _questions = """ 
    test questions of class super
    """

    def __init__(self):
        pass

class test(testsuper):
    _questions = testsuper._questions + "This is the added string"

    def __init__(self):
        print(test._questions)
        pass

if __name__=="__main__":
    t = test()
