from pysurf.qctools import generate_filereader
from pysurf.qctools import Event
from pysurf.qctools import register_event_type
from pysurf.qctools.events import join_events


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
