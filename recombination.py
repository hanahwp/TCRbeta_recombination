import random
from random import choice
from data_arrange import *

base = ("a", "t", "g", "c")
v_insmax = input("Type in the number of maximum insertion for VD junction: ")
j_insmax = input("Type in the number of maximum insertion for DJ junction: ")

def recombinate(TRBLV, TRBD, TRBJ, TRBC):
    #choose 1or2
    n = random.randint(0,1)
    D = TRBD[n] 
    J = TRBJ[n]
    C = TRBC[n]
    V = TRBLV[0]
    
    #choose gene
    vgene = V[choice(V.keys())]
    dgene = D[choice(D.keys())]
    jgene = J[choice(J.keys())]
    cgene = C[choice(C.keys())]

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

    
    #deletion
    vdel_3 = random.randint(0,10)
    jdel_5 = random.randint(0,10)
    ddel_3 = random.randint(0,8)
    ddel_5 = random.randint(0,8)

    #insertion
    vd_insert = "".join([choice(base) for i in  range(random.randint(0,v_insmax))])
    dj_insert = "".join([choice(base) for i in  range(random.randint(0,j_insmax))])

    #recombination
    recomb_assno= vreg[0]+";"+dreg[0]+";"+jreg[0]+";"+creg[0]
    recomb_name= vreg[1]+";"+dreg[1]+";"+jreg[1]+";"+creg[1] 
    recomb_seq= vreg[2][:len(vreg[2])-vdel_3]+vd_insert+dreg[2][ddel_5:-ddel_3]+dj_insert+jreg[2][jdel_5:]+creg[2]


        
    
    # translate
    aa = recomb_seq.translate()
    if aa[0] != "M":
            print recomb_name, "\n", aa, "\n", recomb_seq
            print "\nV before deletion:\n", vreg[2]
            print "deleted amount:", vdel_3
            print "V:\n", vreg[2][:len(vreg[2])-vdel_3]
            print "\nINS VD:", vd_insert
            print "D:", dreg[2][ddel_5:-ddel_3]
            print "INS DJ:", dj_insert
            print "JC:", jreg[2][jdel_5:]+creg[2], "\n\n"

    return recomb_assno, recomb_name, recomb_seq, aa


                                    
def checkseq((recomb_assno, recomb_name, recomb_seq, aa), tcrb, nsense, frameout):
    #assume that shifted frame(not using the conventional stop codon) is degraded
    if '*' in aa:

        keynamevar_dict("%s|%s"%(recomb_assno, recomb_name),str(recomb_seq), nsense )
        return "nsense"
    else:
        if (len(recomb_seq) % 3) == 0:
            keynamevar_dict("%s|%s"%(recomb_assno, recomb_name),str(recomb_seq),tcrb )
            return "tcrb"
        else:
            keynamevar_dict("%s|%s"%(recomb_assno, recomb_name),str(recomb_seq),frameout )
            return "frameout"
            

