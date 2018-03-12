import random
from random import choice
from data_arrange import *
from Bio import SeqIO
import time
import numpy  as np
from  data_arrange import *
from collections import defaultdict as ddict
import matplotlib.pyplot as plt
import csv
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def frequency_from_read(filename):
    text_file = open(filename, "r")
    lines = text_file.read().split('\t')
    text_file.close()

    number_dict = ddict(int)
    for gene_no in lines:
        number_dict[gene_no] += 1
        
    sig = sum(number_dict.values())
    frequency_dict = dict()
    
    for key in number_dict.keys():
        frequency_dict[int(key)] = float(number_dict[key])/float(sig)
    y = {}
    for key in sorted(frequency_dict.keys(), key=int):
        y[key] = frequency_dict[key]

    probs= y.values()
    keys = y.keys()
    return y


d_insvd = frequency_from_read('rawdata/insVD.txt')
d_insdj = frequency_from_read('rawdata/insDJ.txt')

base = ("a", "t", "g", "c")

def recomb_n_check(D_TRBLV, D_TRBD,  D_TRBJ, D_TRBC,D_cys2_index, D_jphe_index):
    t0= time.time()
    insvd= np.random.choice(d_insvd.keys(), replace=True, p=d_insvd.values())
    insdj= np.random.choice(d_insdj.keys(), replace=True, p=d_insdj.values())
    vdel = 12
    jdel = 12

    t1= time.time()

    #choose btw TRBx1 and TRBx2
    n = random.randint(1,2)
    D = D_TRBD['%d'%n] 
    J = D_TRBJ['%d'%n]
    C = D_TRBC['%d'%n]
    V = D_TRBLV

    #choose a gene
    v_genesegment_name=choice(V.keys()) 
    d_genesegment_name=choice(D.keys())
    j_genesegment_name=choice(J.keys())
    c_genesegment_name=choice(C.keys())
    
    vgene = V[v_genesegment_name] 
    dgene = D[d_genesegment_name]
    jgene = J[j_genesegment_name]
    cgene = C[c_genesegment_name]

    #choose an allele
    v_allele_accessno=choice(vgene.keys())
    d_allele_accessno=choice(dgene.keys())
    j_allele_accessno=choice(jgene.keys())
    c_allele_accessno=choice(cgene.keys())

    T_vall=  vgene[v_allele_accessno] #=(alleleno, functionality, position,  length, seq)
    T_dall=  dgene[d_allele_accessno]
    T_jall=  jgene[j_allele_accessno]
    T_call=  cgene[c_allele_accessno]

    v_allele_name = v_genesegment_name+'*'+T_vall[0]
    d_allele_name = d_genesegment_name+'*'+T_dall[0]
    j_allele_name = j_genesegment_name+'*'+T_jall[0]
    c_allele_name = c_genesegment_name+'*'+T_call[0]
    
    
    #number of nucleotide deletion
    vdel_3 = random.randint(0,vdel)
    jdel_5 = random.randint(0,jdel)

    ddel_5 = random.randint(0,12)
    ddel_3 = random.randint(0,12-ddel_5)

    #make nucleotide insertion sequence
    vd_insert = "".join([choice(base) for i in  xrange(insvd)])
    dj_insert = "".join([choice(base) for i in  xrange(insdj)])

    recomb=dict()
    #recombination
    recomb["assno"]= ";".join([
        v_allele_accessno,d_allele_accessno,
        j_allele_accessno,c_allele_accessno
        ])
    recomb["name"]= ';'.join([
        v_allele_name,d_allele_name,
        j_allele_name,c_allele_name
        ])
    
    recomb["seq"]= ''.join([
        T_vall[4][:len(T_vall[4])-vdel_3],vd_insert,
        T_dall[4][ddel_5:-ddel_3],dj_insert,
        T_jall[4][jdel_5:],T_call[4]
        ])
    recomb['seq'] = Seq(recomb['seq'], generic_dna)
    deletion="-".join(["",str(vdel_3), str(ddel_5),str(ddel_3), str(jdel_5)])
    insertion = "+".join(["",str(len(vd_insert)), str(len(dj_insert))])
    modified = str(-vdel_3-ddel_5-ddel_3-jdel_5+len(vd_insert)+len(dj_insert))
    indel = "".join([insertion, deletion])
    indel = "=".join([indel, modified])

    
    #find V-2nd Cys index
    D_cys2_assnos = D_cys2_index[v_genesegment_name]
    D_jphe_assnos = D_jphe_index[j_genesegment_name]

    if (D_cys2_assnos) and (D_jphe_assnos):
        try:
            cys2_index = int(D_cys2_index[v_genesegment_name][v_allele_accessno])
        except KeyError:
            return recomb,  'unconserved', indel
        try:
            jphe_index = int(D_jphe_index[j_genesegment_name][j_allele_accessno])
        except KeyError:
            return recomb,  'unconserved', indel
        if cys2_index+3+1<=len(T_vall[4][:len(T_vall[4])-vdel_3]) and jphe_index+1>=jdel_5:
            aa = recomb["seq"].translate()  #if V-Cys or J-Phe exists, translate
            if '*' in aa:
                return recomb, "nsense",  indel
            else:
                #if V-Cys exists, it means V has been translated in frame
                if (len(recomb["seq"]) % 3) == 0:
                    m = Seq(T_jall[4][jphe_index:], generic_dna) .translate()
                    n  = len(T_call[4])+len(T_jall[4][jphe_index:])
                    if m in aa[-int(round(n/3))-4:-int(round(len(T_call[4])/3))+4]:
                        return recomb, "tcrb", indel
                    else:
                        return recomb, "frameout", indel   #if J-phe is not inframe, it is unexpressed(frameout for now)
                else:
                    return recomb, "frameout", indel    #assume that shifted frame(not using the usual stop codon) is degraded
        else:
            return recomb, "unconserved", indel   #if V-Cys and J-Phe is not in the sequence, consider it as unexpressed(unconserved)
    else:
        return recomb, "unconserved" , indel


    
