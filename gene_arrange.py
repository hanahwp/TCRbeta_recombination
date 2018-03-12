import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
from data_arrange import *


def make_fpa_fasta(genenamestr, genenamelist):
    PSDO= []
    FUNC= []
    if genenamelist:
        for allele in genenamelist:
            if any(c in allele.description for c in ('|F|','|(F)|','|[F]|', '|ORF|','|(ORF)|','|[ORF]|')) :
                FUNC.append(allele)
            if any(c in allele.description for c in ('|P|', '|(P)|','|[P]|')) :
                PSDO.append(allele)
        if "TRBC" in genenamestr:
            cdict = dict()
            for exons in FUNC:
                all_name = exons.id[exons.id.find("|")+1:nth_occur(exons.id, "|", 2)]
                rdict(all_name,exons, cdict)
            FUNC = []
            exn_order =dict()
            for n in xrange(len(cdict.values())):
                for seqs in cdict.values()[n]: #arrange exons according to their order
                    i = int(seqs.description[nth_occur(seqs.description,"|", 5)-1])
                    exn_order[i] = seqs
                    #make joined sequence and new description
                    ntseq = []
                    ntpos=[]
                    ntnt = []
                    for j in sorted(exn_order.keys()):
                        ntseq.append(str(exn_order[j].seq))
                        ntpos.append(exn_order[j].description[nth_occur(exn_order[j].description,"|", 5)+1:nth_occur(exn_order[j].description,"|", 6)])
                        ntnt.append(exn_order[j].description[nth_occur(exn_order[j].description,"|", 6)+1:nth_occur(exn_order[j].description,"|", 7)-2])

                    seq = "".join(ntseq)[1:]
                    position = ";".join(ntpos)
                    region = "C-REGION"
                    nt = "+".join(ntnt).replace(" ", "")

                fle = SeqRecord(Seq(seq, SingleLetterAlphabet()),
                      id=cdict.values()[n][0].id, name=cdict.values()[n][0].name,
                      description="|".join([cdict.values()[n][0].description[:nth_occur(cdict.values()[n][0].description, "|",4)],
                                            region,
                                            position,
                                            nt,
                                            " | | | | | | | |"
                                            ]), 
                      dbxrefs=cdict.values()[n][0].dbxrefs)
                FUNC.append(fle)
                
        SeqIO.write(FUNC, "extdata/%s_F.fasta"%genenamestr, "fasta")
        SeqIO.write(PSDO, "extdata/%s_P.fasta"%genenamestr, "fasta")
        SeqIO.write(genenamelist, "extdata/%s.fasta"%genenamestr, "fasta")  
        
def extract_data(rdir, LVEXON):
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
    TRBD=[]
    TRBJ=[]
    TRBC=[]
    TRBL=[]
    
    for element in TRB:
        if 'TRBV' in element.id:
            if 'L-PART1+L-PART2' in element.description:
                TRBL.append(element)
            else:
                TRBV.append(element)
        elif 'TRBD' in element.id:
            TRBD.append(element)    
        elif 'TRBJ' in element.id:
            TRBJ.append(element)
        elif 'TRBC' in element.id:
            TRBC.append(element)
    make_fpa_fasta("TRBV", TRBV)               
    make_fpa_fasta("TRBD", TRBD)
    make_fpa_fasta("TRBJ", TRBJ)
    make_fpa_fasta("TRBC", TRBC)
    make_fpa_fasta("TRBL", TRBL)
    
    lead = list(SeqIO.parse(LVEXON, "fasta"))
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
    make_fpa_fasta("TRBLV", TRBLV)

