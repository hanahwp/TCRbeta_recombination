import random
from random import choice
from data_arrange import *
from gene_arrange import find_cysindex, find_pheindex

base = ("a", "t", "g", "c")


def recomb_n_check(TRBLV, TRBD, TRBJ, TRBC, v_insmax, j_insmax, vdel, jdel, cys2, jphe):
    #choose btw TRBx1 and TRBx2
    n = random.randint(0,1)
    D = TRBD[n] 
    J = TRBJ[n]
    C = TRBC[n]
    V = TRBLV[0]
    
    #choose a gene
    vgene = V[choice(V.keys())]
    dgene = D[choice(D.keys())]
    jgene = J[choice(J.keys())]
    cgene = C[choice(C.keys())]

    #choose an allele
    call=  choice(cgene)
    dall = choice(dgene)
    jall = choice(jgene)
    vall = choice(vgene)

    #get accesssion no., allele name, and seq
    creg  = {
        "access": call.id[:call.id.find("|")],
        "name" : call.id[call.id.find('TRB'):call.id.find('*')+3],
        "seq": call.seq
        }
    dreg  = {
        "access": dall.id[:call.id.find("|")],
        "name": dall.id[dall.id.find('TRB'):dall.id.find('*')+3],
        "seq": dall.seq
        }
    jreg  = {
        "access":jall.id[:call.id.find("|")],
        "name": jall.id[jall.id.find('TRB'):jall.id.find('*')+3],
        "seq":jall.seq
        }
    vreg  = {
        "access":vall.id[:call.id.find("|")],
        "name":vall.id[vall.id.find('TRB'):vall.id.find('*')+3],
        "seq":vall.seq
        }

    #number of nucleotide deletion
    vdel_3 = random.randint(0,vdel)
    jdel_5 = random.randint(0,jdel)
    ddel_3 = random.randint(0,len(dreg["seq"])/2)
    ddel_5 = random.randint(0,len(dreg["seq"])/2)
    
    #make nucleotide insertion sequence
    vd_insert = "".join([choice(base) for i in  range(random.randint(0,v_insmax))])
    dj_insert = "".join([choice(base) for i in  range(random.randint(0,j_insmax))])

    recomb=dict()
    #recombination
    recomb["assno"]= vreg["access"]+";"+dreg["access"]+";"+jreg["access"]+";"+creg["access"]
    recomb["name"]= vreg["name"]+";"+dreg["name"]+";"+jreg["name"]+";"+creg["name"] 
    recomb["seq"]= vreg["seq"][:len(vreg["seq"])-vdel_3]+vd_insert+dreg["seq"][ddel_5:-ddel_3]+dj_insert+jreg["seq"][jdel_5:]+creg["seq"]
  
    #find V-2nd Cys index
    cysindex = find_cysindex(vall, cys2)
    pheindex = find_pheindex(jall, jphe)
    #if V-Cys and J-Phe is not in the sequence, consider it as unexpressed(unconserved)
    if (cysindex and pheindex) != None:
        if cysindex+3+1<=len(vall.seq[:len(vall.seq)-vdel_3]) and pheindex+1>=jdel_5:
            aa = recomb["seq"].translate()  #if V-Cys or J-Phe exists, translate            
            if '*' in aa:
                return recomb, "nsense"
            else:
                #if V-Cys exists, it means V has been translated inframe
                if (len(recomb["seq"]) % 3) == 0:
                    m =  jreg["seq"][pheindex:].translate()
                    n  = len(call.seq)+len(jreg["seq"][pheindex:])
#                    print "\n",m, "\n",aa[-int(round(n/3))-4:-int(round(len(call.seq)/3))+4], "\n",vall.id,dall.id,jall.id, call.id,"\n", aa
                    if m in aa[-int(round(n/3))-4:-int(round(len(call.seq)/3))+4]:
                        return recomb, "tcrb"
                        print "4"
                    else:
                        return recomb, "frameout"   #if J-phe is not inframe, it is unexpressed(frameout for now)
                        print "5"
                else:
                    return recomb, "frameout"    #assume that shifted frame(not using the usual stop codon) is degraded
        else:
            return recomb, "unconserved"
    else:
        return recomb, "unconserved" 

