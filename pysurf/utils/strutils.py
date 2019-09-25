import re

def split_str(a):
    a = a.strip('(){}[] ')
    asp = re.split(r'[;|,\s\(\)|\*|\n]',a)
    asp_final = []
    for i in range(len(asp)):
        if asp[i] != '':
            asp_final += [asp[i]]
    return asp_final
