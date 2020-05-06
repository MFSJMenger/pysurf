from itertools import permutations
import numpy as np
from copy import deepcopy

class NGridIterator():
    """ Iterator that runs through an n-dimensional grid on spheres around the
        origin, given vectors for the grid points. Thus it runs through an N-
        dimensional grid staying as close as possible to the origin, without 
        visiting points twice.
    """

    def __init__(self, dim):
        """
        Args:
            dim, int
                dim is the dimension of the grid
        """

        self.dim = dim
        # self.R2 is the current square radius of the sphere
        self.R2 = 0
        # self.vectors contains all the vectors on the 
        # sphere with radius sqrt(R2)
        self.vectors = get_vectors(self.dim, self.R2)
        self.vectors_idx = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.vectors_idx >= len(self.vectors):
            self.vectors_idx = 0
            self.vectors = []
            while len(self.vectors) == 0:
                self.R2 += 1
                self.vectors = get_vectors(self.dim, self.R2)
            
        
        res = self.vectors[self.vectors_idx]
        self.vectors_idx += 1

        return res
        



def get_base_vector(dim, maximum, numbers=None):
    """ base vector is a recursive function that 
        calculates a possible combination out of a set,
        given in numbers, for a specific maximum number.
        For example: maximum = 5 and numbers = [4, 1, 1, 0, 0]
        and dim = 2
        Then it will return the vector [4, 1], such that the sum
        is maximum. It will only use numbers which are given in numbers

    """
    if maximum == 0 and dim == 0: return []
    if dim == 0: return False
    if numbers is None: numbers = get_set(dim, maximum)
    vector = []
    
    nb = max(numbers)
    numbers.remove(nb)
    if nb <= maximum:
        vector += [nb]
        res = get_base_vector(dim-1, maximum-nb, numbers)
        if res is False:
            return False
        else:
            vector += res
            return vector
    else:
        res = get_base_vector(dim, maximum, numbers)
        if res is False:
            return False
        else:
            vector += res
            return vector

def get_set(dim, maximum):
    """ Provides all possible square numbers that could be needed
    to calculate maximum as a sum of a dim-dimensional vector.
    For example: dim = 2, maximum = 5, then the function will return
    [4, 1, 1, 0, 0]
    the repetitions of the numbers are important for get_base_vector
    """

    i = 0
    numbers = []
    while i**2 <= maximum:
        n = i**2
        counter = 0
        while n <= maximum and counter < dim:
            numbers += [i**2]
            n += i**2
            counter += 1
        i += 1
    return numbers

def mysum(vec):
    if len(vec) == 0:
        return 0
    else:
        return sum(vec)

def get_vectors(dim, R2):
    """function collecting all vectors for a specific radius
    """

    #collecting base vectors
    base_vecs = []
    numbers = get_set(dim, R2)
    while len(numbers) >= dim:
        vec = get_base_vector(dim, R2, deepcopy(numbers))
        if vec is not False:
            base_vecs += [np.sqrt(vec)]
        numbers.remove(max(numbers))
    #permuting base vectors
    uvecs = []
    for vec in base_vecs:
        for per_vec in permutations(vec):
            uvecs += [per_vec]
    uvecs = list(set(uvecs))

    #adding all possible sign options
    vecs = []
    for vec in uvecs:
        for sign in sign_possibilities(dim):
            vecs += [tuple([int(a*b) for a, b in zip(sign, vec)])]
    vecs = list(set(vecs))
    return vecs


def sign_possibilities(dim):
    """
    giving all possible sign combinations for a specific dimension
    e.g. for dim = 2:
    [[1,1], [-1,1], [1,-1], [-1,-1]]
    """
    vecs = []
    for i in range(dim+1):
        vec = np.ones(dim)
        vec[:i] *= -1
        for svec in permutations(vec):
            vecs += [svec]
    return list(set(vecs))


if __name__=="__main__":
    mygrid = NGridIterator(3)
    for i in range(100):
        print(next(mygrid))
