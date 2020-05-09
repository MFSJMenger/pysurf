from collections import UserDict
import numpy as np


class RequestGenerator:
    """Abstraction to generate Requests consistently"""

    def __init__(self, nstates, properties=None, use_db=False):
        self.nstates = nstates
        # properties that are always asked for (database)
        if properties is None:
            properties = []
        #
        self._request_always = properties
        #
        if use_db is True:
            self.request = self._request_all
        else:
            self.request = self._request

    def _request(self, crd, properties, states=None):
        """add more sanity checks!"""
        properties = properties + self._request_always
        if states is None:
            return Request(crd, properties, list(range(self.nstates)))
        return Request(crd, properties, states)

    def _request_all(self, crd, properties, states=None):
        properties = properties + self._request_always
        return Request(crd, properties, list(range(self.nstates)))

class StateData:

    def __init__(self, states, size):
        self._states = states
        self.data = np.empty((len(states), size), dtype=np.double)

    def __setitem__(self, istate, value):
        self.data[istate] = value

    def __getitem__(self, istate):
        return self.data[istate]


class Request(UserDict):

    def __init__(self, crd, properties, states):
        crd = np.array(crd)
        UserDict.__init__(self)
        #
        data = {prop: None for prop in properties}
        # add state loop
        data['states'] = states
        # add crd
        data['crd'] = crd
        if 'gradient' in data:
            data['gradient'] = StateData(states, crd.size)
        self.data.update(data)
        # store properties
        self.properties = properties
