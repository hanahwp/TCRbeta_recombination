
from collections import Counter
from collections import defaultdict as ddict
import os

def makerdir(insidetcr):
    directory = [os.getcwd(), insidetcr]
    final_directory="".join(directory)
    return final_directory


def rdict(key,value,aDict):
    if not key in aDict:
        aDict[key] = [(value)]
    else:
        aDict[key].append((value))


def keynamevar_dict(key, value, dictionary):
    int = 1
    while int >0:
        if "%s;;;%d" %(key,int) not in dictionary:
            dictionary["%s;;;%d" % (key, int)] = [(value)]
            int= 0
        else:
            int +=1

def most_common(lst):
    data = Counter(lst)
    return max(lst, key=data.get)
