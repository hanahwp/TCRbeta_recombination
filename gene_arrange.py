"""extract VDJC"""
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from data_arrange import *

def extract_data(rdir):
    final_directory = makerdir('/extdata')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)
    rawdata = list(SeqIO.parse(rdir, "fasta"))
    
    homo=[]
    for element in rawdata:
        if 'Homo' in element.id:
            homo.append(element)    
    TRB =[]
    for element in homo:
        if 'TRB' in element.id:
            TRB.append(element)
    SeqIO.write(TRB, "extdata/homosapiense_trb.fasta", "fasta")
    TRBV=[]
    TRBD1=[]
    TRBD2=[]
    TRBJ1=[]
    TRBJ2=[]
    TRBC1=[]
    TRBC2=[]
    TRBV_F = [] 
    TRBV_P = []
    TRBD1_F = []
    TRBD1_P = []
    TRBD2_F = []
    TRBD2_P = []
    TRBJ1_F = []
    TRBJ1_P = []
    TRBJ2_F = []
    TRBJ2_P = []
    TRBC1_F = []
    TRBC1_P = []
    TRBC2_F = []
    TRBC2_P = []
    TRBL = []
    TRBL_F = []
    TRBL_P = []
    
    for element in TRB:
        if 'TRBV' in element.id:
            if 'L-PART1+L-PART2' in element.description:
                TRBL.append(element)
            else:
                TRBV.append(element)
        elif 'TRBD1' in element.id:
            TRBD1.append(element)    
        elif 'TRBD2' in element.id:
            TRBD2.append(element)
        elif 'TRBJ1' in element.id:
            TRBJ1.append(element)
        elif 'TRBJ2' in element.id:
            TRBJ2.append(element)
        elif 'TRBC1' in element.id:
            TRBC1.append(element)
        elif 'TRBC2' in element.id:
            TRBC2.append(element)            
                
    if TRBV:
        for allele in TRBV:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBV_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBV_P.append(allele)
        SeqIO.write(TRBV_F, "extdata/TRBV_F.fasta", "fasta")
    if TRBD1:
        for allele in TRBD1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBD1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBD1_P.append(allele) 
        SeqIO.write(TRBD1_F, "extdata/TRBD1_F.fasta", "fasta")        
    if TRBD2:
        for allele in TRBD2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBD2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBD2_P.append(allele)
        SeqIO.write(TRBD2_F, "extdata/TRBD2_F.fasta", "fasta")        
    if TRBJ1:
        for allele in TRBJ1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBJ1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBJ1_P.append(allele)
        SeqIO.write(TRBJ1_F, "extdata/TRBJ1_F.fasta", "fasta")
    if TRBJ2:
        for allele in TRBJ2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBJ2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBJ2_P.append(allele)        
        SeqIO.write(TRBJ2_F, "extdata/TRBJ2_F.fasta", "fasta")       
    if TRBC1:
        for allele in TRBC1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBC1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBC1_P.append(allele)
        f1= SeqRecord(Seq(str(TRBC1_F[0].seq+TRBC1_F[1].seq+TRBC1_F[2].seq+TRBC1_F[3].seq), SingleLetterAlphabet()),id=TRBC1_F[0].id, name=TRBC1_F[0].name, dbxrefs=TRBC1_F[0].dbxrefs)
        f2= SeqRecord(Seq(str(TRBC1_F[4].seq+TRBC1_F[5].seq+TRBC1_F[6].seq+TRBC1_F[7].seq), SingleLetterAlphabet()),id=TRBC1_F[4].id, name=TRBC1_F[4].name, dbxrefs=TRBC1_F[4].dbxrefs)
        TRBC1_F =[f1,f2]
        SeqIO.write(TRBC1_F, "extdata/TRBC1_F.fasta", "fasta")       
    if TRBC2:
        for allele in TRBC2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBC2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBC2_P.append(allele)
        f1= SeqRecord(Seq(str(TRBC2_F[0].seq+TRBC2_F[1].seq+TRBC2_F[2].seq+TRBC2_F[3].seq), SingleLetterAlphabet()),id=TRBC2_F[0].id, name=TRBC2_F[0].name, dbxrefs=TRBC2_F[0].dbxrefs)
        f2= SeqRecord(Seq(str(TRBC2_F[4].seq+TRBC2_F[5].seq+TRBC2_F[6].seq+TRBC2_F[7].seq), SingleLetterAlphabet()),id=TRBC2_F[4].id, name=TRBC2_F[4].name, dbxrefs=TRBC2_F[4].dbxrefs)
        TRBC2_F =[f1,f2]
        SeqIO.write(TRBC2_F, "extdata/TRBC2_F.fasta", "fasta")       
    if TRBL:
        for allele in TRBL:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBL_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBL_P.append(allele)
        SeqIO.write(TRBL_F, "extdata/TRBL_F.fasta", "fasta")

          
#    print "V allele: ",len(TRBV), "\nD allele: ", len(TRBD1), "+", len(TRBD2),"\nJ allele: ", len(TRBJ1), "+",len(TRBJ2)



    #print "list of [VF,VP], [D1F,D1P], [D2F,D2P],[J1F,J1P], [J2F,J2P],  [C1F,C1P], [C2F,C2P] returned"
    #return [TRBV_F,TRBV_P], [TRBD1_F,TRBD1_P], [TRBD2_F,TRBD2_P], [TRBJ1_F,TRBJ1_P], [TRBJ2_F,TRBJ2_P],[TRBC1_F,TRBC1_P], [TRBC2_F,TRBC2_P]
#    print "C allele: ", len(TRBC1_F), "+",len(TRBC2_F)
              



def l_match(rdir):
    lead = list(SeqIO.parse(rdir, "fasta"))
    TRBV_F = list(SeqIO.parse("extdata/TRBV_F.fasta", "fasta"))
    TRBV_P = list(SeqIO.parse("extdata/TRBV_P.fasta", "fasta"))
    TRBV = TRBV_F+TRBV_P
#    print "total allele: ",  len(TRBV), "\ntype of gene: ", len(set([item.id[item.id.find("TRB"):item.id.find("*")] for item in TRBV]))
    matchedV = []
    match = dict()
    TRBLV =[]
    TRBLV_F = []
    TRBLV_P = []
    
    nostartcodon= 0
    yesmatch = 0
    seqbtwLV = 0
    for allele in TRBV:
        for lvsq in lead:
            if allele.id[allele.id.find("TRB"):allele.id.find("*")+3] in lvsq.id:
                yesmatch +=1
                if allele.seq in lvsq.seq[-len(allele.seq):]:
                    if str(lvsq.seq)[:3] == "atg":
                        temp = [str(lvsq.seq[:-len(allele.seq)]), str(allele.seq)]
                        TRBLV.append(SeqRecord(Seq(str("".join(temp)), SingleLetterAlphabet()), description=allele.description,id=allele.id, name=allele.name, dbxrefs=allele.dbxrefs))
                        rdict(allele.id[allele.id.find("TRB"):allele.id.find("*")],str(lvsq.seq[:-len(allele.seq)]), match)
                        matchedV.append(allele)
                    else: nostartcodon +=1
                else: seqbtwLV +=1
    mismatch= set(TRBV)-set(matchedV)
#    print "total: ", len(set(TRBV)), "found match:", yesmatch, "matchedwithstartcodon:", len(set(matchedV)), "match but no start codon: ", nostartcodon

    #use other allele's lead sequence(are quite similar)
    for allele in mismatch:
        genename = allele.id[allele.id.find("TRB"):allele.id.find("*")]
        for key in match:
            if genename == key:
                temp = [str(most_common(match[key])),str(allele.seq)]
                TRBLV.append(SeqRecord(Seq("".join(temp), SingleLetterAlphabet()), description=allele.description,id=allele.id, name=allele.name, dbxrefs=allele.dbxrefs))
                matchedV.append(allele)
#    print "match:", len(set(matchedV)), "total: ", len(set(TRBV))
    mismatch= set(TRBV)-set(matchedV)
    n = []
    for element in mismatch:
        n.append(element.description[:element.description.find("*")+3])
    print "matching L not found:", len(n)
    print ";".join(n)
    for element in TRBLV:
        if any(c in element.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
            TRBLV_F.append(element)
        if any(c in element.description for c in ('|P|', '|(P)|','|[P]|')) :
            TRBLV_P.append(element)

    SeqIO.write(TRBLV_F, "extdata/TRBLV_F.fasta", "fasta")
    SeqIO.write(TRBLV_P, "extdata/TRBLV_P.fasta", "fasta")


def genedict(inputlist): 
#genedict(TRBJ1_F) -> {'TRBJ1-1':[*1,*2], 'TRBJ1-2':[*1,*2, *3]}
    a = {}
    for j in inputlist:
        a1 =[]
        name = j.id[j.id.find('TRB'):j.id.find('*')]
        for i in inputlist:
            match = i.id[i.id.find('TRB'):i.id.find('*')]
            if name == match:
                a1.append(i)
                a[name] = a1
    yield a


