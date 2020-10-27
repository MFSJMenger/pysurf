from collections.abc import Mapping
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

    def _request(self, crd, properties, states=None, same_crd=False):
        """add more sanity checks!"""
        properties = properties + self._request_always
        if states is None:
            return Request(crd, properties, list(range(self.nstates)), same_crd=same_crd)
        return Request(crd, properties, states, same_crd=same_crd)

    def _request_all(self, crd, properties, states=None, same_crd=False):
        properties = properties + self._request_always
        return Request(crd, properties, list(range(self.nstates)), same_crd=same_crd)


class StateData:

    def __init__(self, states, shape):
        self._states = states
        sh = tuple([len(states)] + list(shape))
        self.data = np.empty(sh, dtype=np.double)

    def set_data(self, data):
        """try to set everything"""
        data = np.asarray(data)
        data = data.reshape(self.data.shape)
        self.data[:] = data

    def __setitem__(self, istate, value):
        idx = self._states.index(istate)
        self.data[idx] = value

    def __getitem__(self, istate):
        idx = self._states.index(istate)
        return self.data[idx]


class Request(Mapping):

    def __init__(self, crd, properties, states, same_crd=False):
        self._properties = {prop: None for prop in properties if prop != 'crd'}
        self.states = states
        self.crd = np.array(crd)
        self.same_crd = same_crd
        #
        if 'gradient' in properties:
            self._properties['gradient'] = StateData(states, self.crd.shape)

    def set(self, name, value):
        """Ignore properties that are not requested!"""
        if name not in self._properties:
            return
        prop = self._properties[name]
        if isinstance(prop, StateData):
            self._set_state_dictionary(prop, value)
        else:
            self._properties[name] = value

    def __getitem__(self, key):
        return self._properties[key]

    def __len__(self):
        return len(self._properties)

    def __iter__(self):
        return iter(self._properties)

    def iter_data(self):
        """Iterate over all data in the request dct"""
        for key, value in self._properties.items():
            if isinstance(value, StateData):
                yield key, value.data
            else:
                yield key, value

    def _set_state_dictionary(self, prop, dct):
        """Set stateData"""
        if not isinstance(dct, Mapping):
            prop.set_data(dct)
            return
        #
        for state, value in dct.items():
            try:
                prop[state] = value
            except ValueError:
                pass
