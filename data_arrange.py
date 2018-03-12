
from collections import Counter
import numpy
import os
import re
#make relative directory
def makerdir(insidetcr):
    directory = [os.getcwd(), insidetcr]
    final_directory="".join(directory)
    return final_directory

#if key exists, add input value
def rdict(key,value,diction):
    if not key in diction:
        diction[key] = [(value)]
    else:
        diction[key].append((value))

            

#number the variation for the same gene combination
def keynamevar_dict(key, value, dictionary):
    int = 1
    while int >0:
        if "%s|%d" %(key,int) not in dictionary:
            dictionary["%s|%d" % (key, int)] = [(value)]
            int= 0
        else:
            int +=1


#count the most common value. If there is a tie, choose the one with the fastest index
def most_common(lst):
    data = Counter(lst)
    return max(lst, key=data.get)

#find the nth occurance for a letter
def nth_occur(string, lookingfor, n, direction=0):
    if direction is 0:
        pos = -1
        for x in xrange(n):
            pos = string.find(lookingfor, pos+1)
            if pos == -1:
                return None
        return pos
    
    if direction is -1:
        pos = len(string)
        for x in xrange(n):
            pos =  string.rfind(lookingfor, 0 ,pos)
            if pos == -1:
                return None
        return pos



def ListToFormattedString(alist):
    # Each item is right-adjusted, width=3
    formatted_list = ['{:>3}' for item in alist]
    s = ''.join(formatted_list)
    return s.format(*alist)
