from qctools import generate_filereader
from qctools import Event
from qctools import register_event_type
from qctools.events import join_events


def event_getter_between():
    
    keyword_args = {}
    args = ['start', 'end']

    def between(iterator, start, end):

        
        if start is None:
            inbetween = True
        else:
            inbetween = False
            # search for start
            for line in iterator:
                if start in line:
                    inbetween = True
                    break
        # return
        if inbetween is False:
            return None, -1
        
        out = []
        for line in iterator:
            if end in line:
                break
            out.append(line)
        return out, 1

    return keyword_args, args, between


register_event_type('between', event_getter_between)


def length(iterator):
    return len(iterator)

NAtoms = Event('NAtoms', 
               'between', {'start': '[Atoms]',
                           'end': '[FREQ]',
                           'ishift': 1},
              func=length,
)

NFreqs = Event('NFreqs', 
               'between', {'start': None,
                           'end': '[FR-COORD]'},
                func=length,
)

Info = join_events(NAtoms, NFreqs)
Info._settings['nmax'] = 1

Freqs = Event('Freqs', 
               'xgrep', {'keyword': '[FREQ]',
                        'ilen': 'NFreqs',
                        'ishift': 1,},
                func='split',
                func_kwargs={'idx': 0, 'typ': float},
                settings = {'reset': True}
)


FrCoords = Event('FrCoords', 
              'xgrep', {'keyword': '[FR-COORD]',
                       'ilen': 'NAtoms',
                       'ishift': 1},
              func='split',
              func_kwargs={
                  'idx': [0, 1, 2, 3], 
                  'typ': [str, float, float, float],
              },
)


def parse_fr_norm_coords(result):
    vibration = {}
    active = -1
    for line in result:
        if 'vibration' in line:
            active += 1
            vibration[active] = []
            continue
        vibration[active].append(list(map(float, line.split())))
    return vibration


FrNormCoords = Event('FrNormCoords', 
        'xgrep', {'keyword': '[FR-NORM-COORD]',
                  'ilen': 'NFreqs*(NAtoms+1)',
                  'ishift': 1,
                  },
        func=parse_fr_norm_coords,
)


MoldenParser = generate_filereader('MoldenParser', {
    'Info': Info,
    'Freqs': Freqs,
    'FrCoords': FrCoords,
    'FrNormCoords': FrNormCoords,
})


class FileIterator:

    def __init__(self, filename):
        with open(filename, 'r') as fh:
            self.lines = fh.readlines()
        self._current = 0
        self._nele = len(self.lines)

    def inc(self):
        self._current += 1

    def next(self):
        if self._current < self._nele:
            line = self.lines[self._current]
            self._current += 1
            return line
        return None
        
    def __iter__(self):
        return self

    def __next__(self):
        line = self.next()
        if line is None:
            raise StopIteration
        return line

    def peek(self, n=0):
        number = self._current + n
        if number < self._nele:
            line = self.lines[number]
            return line

            
def parse_freqs(fh):

    freqs = []

    while True:
        line = fh.peek()
        if line is None or '[' in line: 
            break
        fh.inc()
        freqs.append(float(line))

    return freqs


def parse_vibration(fh):

    vibrations = []
    while True:
        line = fh.peek()
        if line is None or '[' in line or 'vibration' in line: 
            break
        fh.inc()
        x, y, z = line.split()
        vibrations.append([float(x), float(y), float(z)])
    return vibrations


def parse_vibrations(fh):

    vibrations = []
    while True:
        line = fh.peek()
        if 'vibration' in line:
            fh.inc()
            vibrations.append(parse_vibration(fh))
        else:
            break
    return vibrations


def parse_frcoords(fh):
    coords = []
    while True:
        line = fh.peek()
        if '[' in line or line is None:
            break
        atom, x, y, z = line.split()
        coords.append([atom, float(x), float(y), float(z)])
        fh.inc()
    return coords


def parse_molden(filename):

    fh = FileIterator(filename)

    for line in fh:
        line = line.strip()
        if line == '[FREQ]':
            freqs = parse_freqs(fh)
        elif line == '[FR-COORD]':
            coords = parse_frcoords(fh)
        elif line == '[FR-NORM-COORD]':
            vibrations = parse_vibrations(fh)
        #elif '[Atoms]' in line:
        #    coords = parse_coords(line, fh)

    return freqs, coords, vibrations
