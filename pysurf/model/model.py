from abc import ABC, abstractmethod


class Model(ABC):

    @abstractmethod
    def get(self):
        pass
