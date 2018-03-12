
from Bio import SeqIO
from data_arrange import *
from collections import defaultdict as ddict
import csv
import os
from numpy import around
import time


#make a frequency dictionary 
def frequency_dictionary(itemset, total_no_ddict_int, wanted_no_ddict_int):
    #two dictionaries should be defaultdictionary(int)
    outputdict = dict()
    for name in itemset:
        if total_no_ddict_int[name] == 0:
            outputdict[name] = -1 #-100 means not used at all
        else:
            outputdict[name] = numpy.around(float(wanted_no_ddict_int[name])/float(total_no_ddict_int[name]), decimals=5)
    return outputdict



def calculate_percentage(tcrbpath):
    t = int(tcrbpath[8:18])
    tcrb = list(SeqIO.parse(tcrbpath, "fasta"))
    TRBLV = list(SeqIO.parse("extdata/TRBLV_F.fasta", "fasta"))
    TRBD = list(SeqIO.parse("extdata/TRBD_F.fasta", "fasta"))
    TRBJ = list(SeqIO.parse("extdata/TRBJ_F.fasta", "fasta"))
    TRBC = list(SeqIO.parse("extdata/TRBC_F.fasta", "fasta"))
    nsense = open("extdata/nsense.txt").read().splitlines()
    frameout = open("extdata/frameout.txt").read().splitlines()
    unconserved = open("extdata/unconserved.txt").read().splitlines()

    vname=[]
    dname=[]
    jname=[]
    cname=[]

    functionality={}

    for allele in TRBLV:
        vname.append(allele.description.split('|')[1])
        functionality[allele.description.split('|')[1]] = allele.description.split('|')[3]
    for allele in TRBD:
        dname.append(allele.description.split('|')[1])
        functionality[allele.description.split('|')[1]] = allele.description.split('|')[3]
    for allele in TRBJ:
        jname.append(allele.description.split('|')[1])
        functionality[allele.description.split('|')[1]] = allele.description.split('|')[3]
    for allele in TRBC:
        cname.append(allele.description.split('|')[1])
        functionality[allele.description.split('|')[1]] = allele.description.split('|')[3]

    vno_p = ddict(int)
    dno_p = ddict(int)
    jno_p = ddict(int)
    cno_p = ddict(int)
    for recomb in tcrb:
        allname = recomb.description.split('|')[1]
        v,d,j,c = allname.split(';')
        vno_p[v] +=1
        dno_p[d] +=1
        jno_p[j] +=1
        cno_p[c] +=1

    vno_x = ddict(int)
    dno_x = ddict(int)
    jno_x = ddict(int)
    cno_x = ddict(int)
    for item in nsense:
        allname = item.split('|')[1]
        v,d,j,c = allname.split(';')
        vno_x[v] +=1
        dno_x[d] +=1
        jno_x[j] +=1
        cno_x[c] +=1
    for item in frameout:
        allname = item.split('|')[1]
        v,d,j,c = allname.split(';')
        vno_x[v] +=1
        dno_x[d] +=1
        jno_x[j] +=1
        cno_x[c] +=1
    for item in unconserved:
        allname = item.split('|')[1]
        v,d,j,c = allname.split(';')
        vno_x[v] +=1
        dno_x[d] +=1
        jno_x[j] +=1
        cno_x[c] +=1

    vno=ddict(int)
    for name in vname:
        vno[name] =vno_p[name]+vno_x[name] 
    dno=ddict(int)
    for name in dname:
        dno[name] =dno_p[name]+dno_x[name]

    jno=ddict(int)
    for name in jname:
        jno[name] =jno_p[name]+jno_x[name]
    cno=ddict(int)
    for name in cname:
        cno[name] =cno_p[name]+cno_x[name]

    vprdctv = frequency_dictionary(vname, vno, vno_p)
    dprdctv = frequency_dictionary(dname, dno, dno_p)
    jprdctv = frequency_dictionary(jname, jno, jno_p)
    cprdctv = frequency_dictionary(cname, cno, cno_p)



        
    #, sum(vno.values()), sum(vno_p.values())
    if not os.path.exists(makerdir('/extdata/%d_%dallele_productivity_file.csv'%(t,sum(vno_p.values())))):
        allele_productivity_file = open(makerdir('/extdata/%d_%dallele_productivity_file.csv'%(t,sum(vno_p.values()))),'w')

    with open(makerdir('/extdata/%d_%dallele_productivity_file.csv'%(t,sum(vno_p.values()))), 'w') as csvfile:

        fieldnames = ["TYPE", "ALLELE", "FUNCTIONALITY", "CHOSEN FREQUENCY", "RESULT_PRODUCTIVITY", "UNIQUE_PRODUCTIVITY"]
        logwriter = csv.DictWriter(allele_productivity_file, fieldnames=fieldnames)
        logwriter.writeheader()

        for allele in vprdctv.keys():
            logwriter.writerow({"TYPE":"V",
                                "ALLELE":allele,
                                "FUNCTIONALITY": functionality[allele],
                                "CHOSEN FREQUENCY":around(float(vno[allele])/float(sum(vno.values())), decimals=5),
                                "RESULT_PRODUCTIVITY":around(float(vno_p[allele])/float(sum(vno_p.values())), decimals=5),
                                "UNIQUE_PRODUCTIVITY":around(float(vno_p[allele])/float(vno[allele]),decimals=5)})
        for allele in dprdctv.keys():
            logwriter.writerow({"TYPE":"D",
                                "ALLELE":allele,
                                "FUNCTIONALITY": functionality[allele],
                                "CHOSEN FREQUENCY":around(float(dno[allele])/float(sum(dno.values())), decimals=5),
                                "RESULT_PRODUCTIVITY":around(float(dno_p[allele])/float(sum(dno_p.values())), decimals=5),
                                "UNIQUE_PRODUCTIVITY":around(float(dno_p[allele])/float(dno[allele]),decimals=5)})        
        for allele in jprdctv.keys():
            logwriter.writerow({"TYPE":"J",
                                "ALLELE":allele,
                                "FUNCTIONALITY": functionality[allele],
                                "CHOSEN FREQUENCY":around(float(jno[allele])/float(sum(jno.values())), decimals=5),
                                "RESULT_PRODUCTIVITY":around(float(jno_p[allele])/float(sum(jno_p.values())), decimals=5),
                                "UNIQUE_PRODUCTIVITY":around(float(jno_p[allele])/float(jno[allele]),decimals=5)})        
        for allele in cprdctv.keys():
            logwriter.writerow({"TYPE":"C",
                                "ALLELE":allele,
                                "FUNCTIONALITY": functionality[allele],
                                "CHOSEN FREQUENCY":around(float(cno[allele])/float(sum(cno.values())), decimals=5),
                                "RESULT_PRODUCTIVITY":around(float(cno_p[allele])/float(sum(cno_p.values())), decimals=5),
                                "UNIQUE_PRODUCTIVITY":around(float(cno_p[allele])/float(cno[allele]),decimals=5)})        

    print "file made : /extdata/%d_%dallele_productivity_file.csv"%(t,sum(vno_p.values())) 
    allele_productivity_file.close()
    return "/extdata/%d_%dallele_productivity_file.csv"%(t,sum(vno_p.values())) 


