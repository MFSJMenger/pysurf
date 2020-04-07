from collections import UserDict


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


class Request(UserDict):

    def __init__(self, crd, properties, states):
        data = {prop: None for prop in properties}
        # add state loop
        data['states'] = states
        # add crd
        data['crd'] = crd
        UserDict.__init__(self)
        self.data.update(data)
        # store properties
        self.properties = properties
