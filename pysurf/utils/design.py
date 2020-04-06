from itertools import cycle
from collections import namedtuple


RGBColor = namedtuple("RGBColor", ["red", "green", "blue"])
#
BLACK = RGBColor(0,0,0)                                                 
BLUE = RGBColor(63/255, 81/255, 181/255)                                
RED = RGBColor(253/255, 86/255, 33/255)                                 
YELLOW = RGBColor(255/255, 233/255, 75/255)                             
GREEN = RGBColor(138/255, 193/255, 73/255)                              
#
colors = cycle([BLACK, BLUE, RED, GREEN, YELLOW])
