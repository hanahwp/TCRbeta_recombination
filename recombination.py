import numpy as np
import random
from random import choice

def recombinate(TRBLV, TRBD, TRBJ, TRBC):
    #choose 1or2
    n = np.random.randint(0,2)
    D = TRBD[n] 
    J = TRBJ[n]
    C = TRBC[n]
    V = TRBLV[0]
    
    #choose gene
    vgene = V[random.choice(V.keys())]
    dgene = D[random.choice(D.keys())]
    jgene = J[random.choice(J.keys())]
    cgene = C[random.choice(C.keys())]

    #choose allele
    call=  choice(cgene)
    dall = choice(dgene)
    jall = choice(jgene)
    vall = choice(vgene)

    #get accesssion no., allele name, and seq
    creg  = [call.id[:6], call.id[call.id.find('TRB'):call.id.find('*')+3], call.seq]
    dreg  = [dall.id[:6], dall.id[dall.id.find('TRB'):dall.id.find('*')+3], dall.seq]
    jreg  = [jall.id[:6], jall.id[jall.id.find('TRB'):jall.id.find('*')+3], jall.seq]
    vreg  = [vall.id[:6], vall.id[vall.id.find('TRB'):vall.id.find('*')+3], vall.seq]

    #recombination
    recomb_assno= vreg[0]+dreg[0]+jreg[0]+creg[0]
    recomb_name= vreg[1]+";"+dreg[1]+";"+jreg[1]+";"+creg[1] 
    recomb_seq= vreg[2]+dreg[2]+jreg[2]+creg[2]
    
    # translate
    aa = recomb_seq.translate()
    return recomb_assno, recomb_name, recomb_seq, aa


def checkseq((recomb_assno, recomb_name, recomb_seq, aa), tcrb, nsense, frameout):
    name = {"tcrb":"tcrb", "nsense":"nsense", "frameout":"frameout"}
    #assume that shifted frame(not using the conventional stop codon) is degraded
    if '*' in aa:
        nsense["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
        return name["nsense"]
    else:
        if (len(recomb_seq) % 3) == 0:
            tcrb["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
            return name["tcrb"]
        else:
            frameout["%s|%s"%(recomb_assno, recomb_name)] = recomb_seq
            return name["frameout"]
            

