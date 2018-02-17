"""extract VDJC"""
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from data_arrange import *

def extract_data(rdir):
    #make /extdata directory
    final_directory = makerdir('/extdata')
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)
        
    #read fasta file
    rawdata = list(SeqIO.parse(rdir, "fasta"))
    homo=[]
    for element in rawdata:
        if 'Homo' in element.id:
            homo.append(element)    
    TRB =[]
    for element in homo:
        if 'TRB' in element.id:
            TRB.append(element)
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
        SeqIO.write(TRBV_P, "extdata/TRBV_P.fasta", "fasta")
        SeqIO.write(TRBV, "extdata/TRBV.fasta", "fasta")
    if TRBD1:
        for allele in TRBD1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBD1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBD1_P.append(allele) 
        SeqIO.write(TRBD1_F, "extdata/TRBD1_F.fasta", "fasta")        
        SeqIO.write(TRBD1_P, "extdata/TRBD1_P.fasta", "fasta")
        SeqIO.write(TRBD1, "extdata/TRBD1.fasta", "fasta")
    if TRBD2:
        for allele in TRBD2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBD2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBD2_P.append(allele)
        SeqIO.write(TRBD2_F, "extdata/TRBD2_F.fasta", "fasta")
        SeqIO.write(TRBD2_P, "extdata/TRBD2_P.fasta", "fasta")
        SeqIO.write(TRBD2, "extdata/TRBD2.fasta", "fasta")
    if TRBJ1:
        for allele in TRBJ1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBJ1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBJ1_P.append(allele)
        SeqIO.write(TRBJ1_F, "extdata/TRBJ1_F.fasta", "fasta")
        SeqIO.write(TRBJ1_P, "extdata/TRBJ1_P.fasta", "fasta")
        SeqIO.write(TRBJ1, "extdata/TRBJ1.fasta", "fasta")
    if TRBJ2:
        for allele in TRBJ2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBJ2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBJ2_P.append(allele)        
        SeqIO.write(TRBJ2_F, "extdata/TRBJ2_F.fasta", "fasta")
        SeqIO.write(TRBJ2_P, "extdata/TRBJ2_P.fasta", "fasta")
        SeqIO.write(TRBJ2, "extdata/TRBJ2.fasta", "fasta")
    if TRBC1:
        for allele in TRBC1:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBC1_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBC1_P.append(allele)
           
        position_1 = "".join([TRBC1_F[0].description[nth_occur(TRBC1_F[0].description, "|", 5)+1:nth_occur(TRBC1_F[0].description, "|", 6)],
                           TRBC1_F[1].description[nth_occur(TRBC1_F[1].description, "|", 5)+2:nth_occur(TRBC1_F[1].description, "|", 6)],
                           TRBC1_F[2].description[nth_occur(TRBC1_F[2].description, "|", 5)+2:nth_occur(TRBC1_F[2].description, "|", 6)],
                           ",",
                           TRBC1_F[3].description[nth_occur(TRBC1_F[3].description, "|", 5)+2:nth_occur(TRBC1_F[3].description, "|", 6)]
                           ])
        reg_name = "C-REGION"
        nt_1 = "+".join([TRBC1_F[0].description[nth_occur(TRBC1_F[0].description, "|", 6)+1:nth_occur(TRBC1_F[0].description, "|", 7)-3],
                       TRBC1_F[1].description[nth_occur(TRBC1_F[1].description, "|", 6)+1:nth_occur(TRBC1_F[1].description, "|", 7)-3],
                       TRBC1_F[2].description[nth_occur(TRBC1_F[2].description, "|", 6)+1:nth_occur(TRBC1_F[2].description, "|", 7)-3],
                       TRBC1_F[3].description[nth_occur(TRBC1_F[3].description, "|", 6)+1:nth_occur(TRBC1_F[3].description, "|", 7)-3],
                       " nt"
                       ])

        position_2 = "".join([TRBC1_F[4].description[nth_occur(TRBC1_F[4].description, "|", 5)+1:nth_occur(TRBC1_F[4].description, "|", 6)],
                           TRBC1_F[5].description[nth_occur(TRBC1_F[5].description, "|", 5)+2:nth_occur(TRBC1_F[5].description, "|", 6)],
                           TRBC1_F[6].description[nth_occur(TRBC1_F[6].description, "|", 5)+2:nth_occur(TRBC1_F[6].description, "|", 6)],
                           ",",
                           TRBC1_F[7].description[nth_occur(TRBC1_F[7].description, "|", 5)+2:nth_occur(TRBC1_F[7].description, "|", 6)]
                           ])
        nt_2 = "+".join([TRBC1_F[4].description[nth_occur(TRBC1_F[4].description, "|", 6)+1:nth_occur(TRBC1_F[4].description, "|", 7)-3],
                       TRBC1_F[5].description[nth_occur(TRBC1_F[5].description, "|", 6)+1:nth_occur(TRBC1_F[5].description, "|", 7)-3],
                       TRBC1_F[6].description[nth_occur(TRBC1_F[6].description, "|", 6)+1:nth_occur(TRBC1_F[6].description, "|", 7)-3],
                       TRBC1_F[7].description[nth_occur(TRBC1_F[7].description, "|", 6)+1:nth_occur(TRBC1_F[7].description, "|", 7)-3],
                       " nt"
                       ])                
            
                
        f1= SeqRecord(Seq(str(TRBC1_F[0].seq[1:]+TRBC1_F[1].seq+TRBC1_F[2].seq+TRBC1_F[3].seq), SingleLetterAlphabet()),
                      id=TRBC1_F[0].id, name=TRBC1_F[0].name,
                      description="|".join([TRBC1_F[0].description[:nth_occur(TRBC1_F[0].description, "|", 4)],
                                           reg_name,
                                           position_1,
                                           nt_1,
                                           TRBC1_F[0].description[nth_occur(TRBC1_F[0].description, "|",7)+1:]
                                           ]), 
                      dbxrefs=TRBC1_F[0].dbxrefs)
        f2= SeqRecord(Seq(str(TRBC1_F[4].seq[1:]+TRBC1_F[5].seq+TRBC1_F[6].seq+TRBC1_F[7].seq), SingleLetterAlphabet()),
                      id=TRBC1_F[4].id, name=TRBC1_F[4].name,
                      description="|".join([TRBC1_F[4].description[:nth_occur(TRBC1_F[4].description, "|", 4)],
                                           reg_name,
                                           position_2,
                                           nt_2,
                                           TRBC1_F[4].description[nth_occur(TRBC1_F[4].description, "|",7)+1:]
                                           ]), 
                      dbxrefs=TRBC1_F[4].dbxrefs)
        TRBC1_F =[f1,f2]
        SeqIO.write(TRBC1_F, "extdata/TRBC1_F.fasta", "fasta")       
    if TRBC2:
        for allele in TRBC2:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBC2_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBC2_P.append(allele)



            
        position_1 = "".join([TRBC2_F[0].description[nth_occur(TRBC2_F[0].description, "|", 5)+1:nth_occur(TRBC2_F[0].description, "|", 6)],
                           TRBC2_F[1].description[nth_occur(TRBC2_F[1].description, "|", 5)+2:nth_occur(TRBC2_F[1].description, "|", 6)],
                           TRBC2_F[2].description[nth_occur(TRBC2_F[2].description, "|", 5)+2:nth_occur(TRBC2_F[2].description, "|", 6)],
                           ",",
                           TRBC2_F[3].description[nth_occur(TRBC2_F[3].description, "|", 5)+2:nth_occur(TRBC2_F[3].description, "|", 6)]
                           ])
        reg_name = "C-REGION"
        nt_1 = "+".join([TRBC2_F[0].description[nth_occur(TRBC2_F[0].description, "|", 6)+1:nth_occur(TRBC2_F[0].description, "|", 7)-3],
                       TRBC2_F[1].description[nth_occur(TRBC2_F[1].description, "|", 6)+1:nth_occur(TRBC2_F[1].description, "|", 7)-3],
                       TRBC2_F[2].description[nth_occur(TRBC2_F[2].description, "|", 6)+1:nth_occur(TRBC2_F[2].description, "|", 7)-3],
                       TRBC2_F[3].description[nth_occur(TRBC2_F[3].description, "|", 6)+1:nth_occur(TRBC2_F[3].description, "|", 7)-3],
                       " nt"
                       ])

        position_2 = "".join([TRBC2_F[4].description[nth_occur(TRBC2_F[4].description, "|", 5)+1:nth_occur(TRBC2_F[4].description, "|", 6)],
                           TRBC2_F[5].description[nth_occur(TRBC2_F[5].description, "|", 5)+2:nth_occur(TRBC2_F[5].description, "|", 6)],
                           TRBC2_F[6].description[nth_occur(TRBC2_F[6].description, "|", 5)+2:nth_occur(TRBC2_F[6].description, "|", 6)],
                           ",",
                           TRBC2_F[7].description[nth_occur(TRBC2_F[7].description, "|", 5)+2:nth_occur(TRBC2_F[7].description, "|", 6)]
                           ])
        nt_2 = "+".join([TRBC2_F[4].description[nth_occur(TRBC2_F[4].description, "|", 6)+1:nth_occur(TRBC2_F[4].description, "|", 7)-3],
                       TRBC2_F[5].description[nth_occur(TRBC2_F[5].description, "|", 6)+1:nth_occur(TRBC2_F[5].description, "|", 7)-3],
                       TRBC2_F[6].description[nth_occur(TRBC2_F[6].description, "|", 6)+1:nth_occur(TRBC2_F[6].description, "|", 7)-3],
                       TRBC2_F[7].description[nth_occur(TRBC2_F[7].description, "|", 6)+1:nth_occur(TRBC2_F[7].description, "|", 7)-3],
                       " nt"
                       ])                
            


        f1= SeqRecord(Seq(str(TRBC2_F[0].seq[1:]+TRBC2_F[1].seq+TRBC2_F[2].seq+TRBC2_F[3].seq), SingleLetterAlphabet()),
                      id=TRBC2_F[0].id, name=TRBC2_F[0].name,
                      description="|".join([TRBC2_F[0].description[:nth_occur(TRBC2_F[0].description, "|", 4)],
                                           reg_name,
                                           position_1,
                                           nt_1,
                                           TRBC2_F[0].description[nth_occur(TRBC2_F[0].description, "|",7)+1:]
                                           ]), 
                      dbxrefs=TRBC2_F[0].dbxrefs)
        f2= SeqRecord(Seq(str(TRBC2_F[4].seq[1:]+TRBC2_F[5].seq+TRBC2_F[6].seq+TRBC2_F[7].seq), SingleLetterAlphabet()),
                      id=TRBC2_F[4].id, name=TRBC2_F[4].name,
                      description="|".join([TRBC2_F[4].description[:nth_occur(TRBC2_F[4].description, "|", 4)],
                                           reg_name,
                                           position_2,
                                           nt_2,
                                           TRBC2_F[4].description[nth_occur(TRBC2_F[4].description, "|",7)+1:]
                                           ]),
                      dbxrefs=TRBC2_F[4].dbxrefs)
        TRBC2_F =[f1,f2]
        SeqIO.write(TRBC2_F, "extdata/TRBC2_F.fasta", "fasta")       
    if TRBL:
        for allele in TRBL:
            if any(c in allele.description for c in ('|F|', '|ORF|','|(F)|','|[F]|','|(ORF)|','|[ORF]|')) :
                TRBL_F.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                TRBL_P.append(allele)
        SeqIO.write(TRBL_F, "extdata/TRBL_F.fasta", "fasta")
        SeqIO.write(TRBL_P, "extdata/TRBL_P.fasta", "fasta")
    print len(TRBV)
    print len(TRBV_P)+len(TRBV_F)



def l_match(rdir):
    lead = list(SeqIO.parse(rdir, "fasta"))
    TRBV_F = list(SeqIO.parse("extdata/TRBV_F.fasta", "fasta"))
    TRBV_P = list(SeqIO.parse("extdata/TRBV_P.fasta", "fasta"))
    TRBV = TRBV_F+TRBV_P

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
                        length ="{0}+{1} nt".format(len(temp[0]), len(temp[1]))
                        new_description = allele.description[:nth_occur(allele.description, "|",6)+1]+length + allele.description[nth_occur(allele.description, "|", 7):]
                        TRBLV.append(SeqRecord(Seq(str("".join(temp)), SingleLetterAlphabet()), description=new_description, id=lvsq.id, name=lvsq.name, dbxrefs=lvsq.dbxrefs))
                        rdict(allele.id[allele.id.find("TRB"):allele.id.find("*")],str(lvsq.seq[:-len(allele.seq)]), match)
                        matchedV.append(allele)
                    else: nostartcodon +=1
                else: seqbtwLV +=1
    mismatch= set(TRBV)-set(matchedV)

    #for V alleles with no matching L, use other allele's L sequence
    for allele in mismatch:
        genename = allele.id[allele.id.find("TRB"):allele.id.find("*")]
        for key in match:
            if genename == key:
                temp = [str(most_common(match[key])),str(allele.seq)]
                length =  "{0}+{1} nt".format(len(temp[0]),len(temp[1]))
                new_description = allele.description[:nth_occur(allele.description, "|",6)+1]+length + allele.description[nth_occur(allele.description, "|", 7):]
                TRBLV.append(SeqRecord(Seq("".join(temp), SingleLetterAlphabet()), description=new_description,id=allele.id, name=allele.name, dbxrefs=allele.dbxrefs))
                matchedV.append(allele)

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
    SeqIO.write(TRBLV,"extdata/TRBLV.fasta","fasta")

    #find 2ndcys index in LV sequence


def find_cysindex(item, cys2):
    for cons in cys2:
        if item.id[:item.id.find("|")] in cons.id:  
            if item.id[item.id.find("RBV"):item.id.find("*")+3] in cons.id:
                vpos =  item.description[nth_occur(item.description, "|", 5)+1:nth_occur(item.description, ".", 2, -1)]
                if 'complement' in cons.description:
                    cyspos = cons.description[cons.description.find("complement")+11:nth_occur(cons.description,".", 2,-1)]
                    cysindex=-(int(cyspos)-int(vpos)+3)
                    return cysindex
                else: 
                    cyspos = cons.description[nth_occur(cons.description, "|", 5)+1:nth_occur(cons.description, ".", 2, -1)]
                    l_length = item.description[nth_occur(item.description, "|",6)+1:nth_occur(item.description, "|",7)]
                    l_length = l_length[:l_length.find("+")]
                    cysindex=int(cyspos)-int(vpos)+int(l_length)
                    return cysindex
            elif item.id[item.id.find("RBV"):item.id.find("*")] in cons.id:
                
                vpos =  item.description[nth_occur(item.description, "|", 5)+1:nth_occur(item.description, ".", 2, -1)]
                if 'complement' in cons.description:
                    cyspos = cons.description[cons.description.find("complement")+11:nth_occur(cons.description,".", 2,-1)]
                    cysindex= -(int(cyspos)-int(vpos)+3)
                    return cysindex
                else:                        
                    cyspos = cons.description[nth_occur(cons.description, "|", 5)+1:nth_occur(cons.description, ".", 2, -1)]
                    l_length = item.description[nth_occur(item.description, "|",6)+1:nth_occur(item.description, "|",7)]
                    l_length = l_length[:l_length.find("+")]
                    cysindex=int(cyspos)-int(vpos)+int(l_length)
                    return cysindex           

def find_pheindex(item, jphe):
    for cons in jphe:
        if item.id[:item.id.find("|")] in cons.id:
            if item.id[item.id.find("RBJ"):item.id.find("*")+3] in cons.id:
                jpos =  item.description[nth_occur(item.description, "|", 5)+1:nth_occur(item.description, ".", 2, -1)]
                if 'complement' in cons.description:
                    phepos = cons.description[cons.description.find("complement")+11:nth_occur(cons.description,".", 2,-1)]
                    pheindex=-(int(phepos)-int(jpos)+3)
                    return pheindex
                else: 
                    phepos = cons.description[nth_occur(cons.description, "|", 5)+1:nth_occur(cons.description, ".", 2, -1)]
                    pheindex=int(phepos)-int(jpos)
                    return pheindex
            elif item.id[item.id.find("RBJ"):item.id.find("*")] in cons.id:
                
                jpos =  item.description[nth_occur(item.description, "|", 5)+1:nth_occur(item.description, ".", 2, -1)]
                if 'complement' in cons.description:
                    phepos = cons.description[cons.description.find("complement")+11:nth_occur(cons.description,".", 2,-1)]
                    pheindex= -(int(phepos)-int(jpos)+3)
                    return pheindex
                else:                        
                    phepos = cons.description[nth_occur(cons.description, "|", 5)+1:nth_occur(cons.description, ".", 2, -1)]
                    pheindex=int(phepos)-int(jpos)
                    return pheindex
                

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




