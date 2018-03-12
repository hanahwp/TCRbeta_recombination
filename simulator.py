
## this program generates translatable rearranged trb genes
## input file is "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.txt", downloaded from IMGT
## output file is [explanation.txt], tcrb.fasta
## 2018. 2. 27.
## by HWPark at CSBL, Korea University

import time
t0=time.time()
import sys
import numpy 
from collections import defaultdict as ddict
from argparse import ArgumentParser
from gene_arrange import *
from recombination import *
from calculate_percentage import *
from k_clustering import *

def main(argv):
    argparse_usage = 'simulator.py -i <input_fasta> -n <wanted_number>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-i", "--input_fasta", dest="input_fasta", nargs=1,
        help="input fasta file"
    )
    parser.add_argument(
        '-n',"--wanted_number", dest='wanted_number', nargs=1,
        help='number of tcr beta that will be generated'
        )
    args = parser.parse_args()
    
    D_input = {}
    if args.input_fasta:
        D_input['input_fasta'] = os.path.abspath(args.input_fasta[0])
    else:
        print '[ERROR] Please provide INPUT FASTA'
        parser.print_help()
        sys.exit(2)
    if args.wanted_number:

        D_input['wanted_number'] = args.wanted_number[0]
    else:
        print '[ERROR] Please type in the number of TCRB you want to generate'
        parser.print_help()
        sys.exit(2)
    return D_input            


if __name__ == "__main__":
    D_input = main(sys.argv[1:])

LVEXON_file = "rawdata/LVEXON.txt"
extract_data(D_input['input_fasta'], LVEXON_file)


def fasta2genedict(file_path):
    with open(file_path, mode='r') as allele_file:
        # FASTA records (thus reducing memory footprint)
        D_gene_allele = ddict(dict)
        for record in SeqIO.parse(allele_file, 'fasta'):
            description = record.description
            sequence = str(record.seq)
            description_split = description.split('|')
            accession_number = description_split[0].replace('>','')
            gene_name = description_split[1].split('*')[0]
            allele_number = description_split[1].split('*')[1]
            functionality =  description_split[3]
            position = description_split[5].split('.')[0]
            length = description_split[6].replace(' nt','')
            D_gene_allele[(gene_name)].update({accession_number:(allele_number,functionality, position, length, sequence)})
            gene_name, accession_number, allele_number,functionality, position, length, sequence = (None,)*7
    return D_gene_allele


TRBLV_file = 'extdata/TRBLV_F.fasta'
D_TRBLV= fasta2genedict(TRBLV_file)


def genecluster_grouping(d_gene):
    D_gene =ddict(dict)
    for key in d_gene.keys():
        D_gene[key[4]].update({key:d_gene[key]})
    return D_gene


TRBD_file = 'extdata/TRBD_F.fasta'
d_TRBD = fasta2genedict(TRBD_file)
D_TRBD = genecluster_grouping(d_TRBD)

TRBJ_file = 'extdata/TRBJ_F.fasta'
d_TRBJ = fasta2genedict(TRBJ_file)
D_TRBJ = genecluster_grouping(d_TRBJ)

TRBC_file = 'extdata/TRBC_F.fasta'
d_TRBC = fasta2genedict(TRBC_file)
D_TRBC = genecluster_grouping(d_TRBC)


# Parse cys2, jphe
def make_cysphe_position_dict(cysphe_file):
    with open(cysphe_file, mode='r') as record_file:
        D_cysphe = ddict(dict)
        for record in SeqIO.parse(record_file, 'fasta'):
            # Extract individual parts of the FASTA record
            description = record.description            
            description_split = description.split('|')
            accession_number = description_split[0].replace('>','').split('.')[0]
            alleles = description_split[1]
            T_alleles = ()
            S_alleles = alleles.replace('or', "*").split('*')
            for alleles_split in S_alleles:
                if alleles_split.startswith('T'):
                    gene_name = alleles_split.replace('TCRBV', 'TRBV')
                else:
                    T_alleles += (alleles_split,)
            position = description_split[5]
            position_split = position.split('..')
            D_cysphe[(gene_name)].update({accession_number: (position_split[0],) + T_alleles})
            gene_name, accession_number, position_split[0], T_alleles =4*(None,)
        return D_cysphe


cys2_file = 'rawdata/cys2.txt'
D_cys2_position = make_cysphe_position_dict(cys2_file)
jphe_file = 'rawdata/jphe.txt'
D_jphe_position = make_cysphe_position_dict(jphe_file)

D_cys2_index = ddict(dict)
for v_genename in D_TRBLV.keys():
    if v_genename in D_cys2_position: #ass_no and genesegment name  is identical
        for accession_number in D_TRBLV[v_genename].keys():
            if accession_number in D_cys2_position[v_genename]:
                v_position = D_TRBLV[v_genename][accession_number][2]
                cys_position = D_cys2_position[v_genename][accession_number][0]
                if 'complement' in cys_position:
                    cys_position=cys_position.split('(')[1]
                    cys_index= -(int(cys_position)-int(v_position)+3)
                    D_cys2_index[v_genename].update({accession_number:cys_index})
                else:
                    l_length = D_TRBLV[v_genename][accession_number][3].split('+')[0]
                    cys_index=int(cys_position)-int(v_position)+int(l_length)
                    D_cys2_index[v_genename].update({accession_number:cys_index})
            else: #thereis no matching ass-noand genesegment-name
                continue #print v_genename
    else:
        continue #print v_genename

D_jphe_index = ddict(dict)
for D_geneclustergroup in D_TRBJ.values():
    for j_genename in D_geneclustergroup.keys():
        if j_genename in D_jphe_position: #ass_no and genesegment name  is identical
            for accession_number in D_geneclustergroup[j_genename].keys():
                if accession_number in D_jphe_position[j_genename]:
                    j_position = D_geneclustergroup[j_genename][accession_number][2]
                    phe_position = D_jphe_position[j_genename][accession_number][0]
                    if 'complement' in phe_position:
                        phe_position=phe_position.split('(')[1]
                        phe_index= -(int(phe_position)-int(j_position)+3)
                        D_jphe_index[j_genename].update({accession_number:phe_index})
                    else:
                        phe_index=int(phe_position)-int(j_position)
                        D_jphe_index[j_genename].update({accession_number:phe_index})
                else: #thereis no matching ass-noand genesegment-name
                    continue #print j_genename
        else:
            continue #print j_genename

tcrb = ddict(list)
nsense = ddict(list)
frameout = ddict(list)
unconserved = ddict(list)
loop_no = 0
tcrb_no = 0
nsense_no = 0
frameout_no = 0
unconserved_no = 0
t1=time.time()
print "preparation time: ",t1-t0

while tcrb_no < int(D_input['wanted_number']):
    t2=time.time()
    #returns no., name, seq, aa
    rearranged = recomb_n_check(D_TRBLV, D_TRBD,D_TRBJ, D_TRBC,D_cys2_index,D_jphe_index)
    t3 = time.time()
#    print "recombination and proof reading:", t3-t2
    if rearranged[1] == "tcrb":
        tcrb["%s|%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"],rearranged[2])].append(str(rearranged[0]["seq"]))
        tcrb_no += 1
    elif rearranged[1] == "frameout":
        frameout["%s|%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"],rearranged[2])].append(str(rearranged[0]["seq"]))
        frameout_no += 1
    elif rearranged[1]  == "nsense":
        nsense["%s|%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"],rearranged[2])].append(str(rearranged[0]["seq"]))
        nsense_no += 1
    elif rearranged[1] == "unconserved":
        unconserved["%s|%s|%s"%(rearranged[0]["assno"], rearranged[0]["name"],rearranged[2])].append(str(rearranged[0]["seq"]))
        unconserved_no+=1
    loop_no +=1
    
t4=time.time()


print "\nTotal:", loop_no#, "\nUnique:", len(tcrb)+len(frameout)+len(nsense)+len(unconserved)
print "Productive_Total:", tcrb_no, len(tcrb)#, "\nProductive_Unique:",  len(tcrb)
print "Out_of_Frame_Total:", frameout_no, len(frameout)#, "\nOut_of_Frame_Unique:", len(frameout)
print "Has_Stop_Total:", nsense_no, len(nsense)#, "\nHas_Stop_Unique:", len(nsense)
print "Unconserved_Total:", unconserved_no, len(unconserved)#, "\nUnconserved_Unique:", len(unconserved)
if loop_no == 0:
    print "Productive_Total/Total % cannot be calculated. (Total = 0) \n"
else:
    print "Productive_Total/Total % : ", numpy.around((100*float(tcrb_no))/float(loop_no),decimals=2), "\n"
with open("extdata/%d_%dtcrb.fasta"%(t0, tcrb_no), "a+") as tcrbfile:
    for combination in tcrb.keys():
        for i in range(len(tcrb[combination])):
            line = [SeqRecord(Seq(tcrb[combination][i],SingleLetterAlphabet()), id=combination+"|%d"%(i+1))]
            SeqIO.write(line, tcrbfile, 'fasta' )

    
nsensetxt = open("extdata/nsense.txt", "w")
for key in nsense:
    nsensetxt.write("%s\n"%key)
nsensetxt.close()
frameouttxt = open("extdata/frameout.txt", "w")
for key in frameout:
    frameouttxt.write("%s\n"%key)
frameouttxt.close()
unconservedtxt = open("extdata/unconserved.txt", "w")
for key in unconserved:
    unconservedtxt.write("%s\n"%key)
unconservedtxt.close()

explanation = open("extdata/explanation_of_trb.fasta.txt", "w")
readme = ("1.'A00000' is IMGT/LIGM-DB accession number. Four accession numbers of each V, D, J, and C genes are serially aligned\n",
          "2. gene and allele names are serially written with semicolon in between\n",
          "3. number after '|' shows the ordinal number of variations generated from the same combination of genes\n",
          "4. combined sequence\n",
          "#by HWPark at CSBL, Korea University")
explanation.write(''.join(readme))
explanation.close()

t5=time.time()
print 'writing time:',t5-t4
print 'total time :', t5-t0


filename = calculate_percentage("extdata/%d_%dtcrb.fasta"%(t0, tcrb_no))
k_clustering(filename[1:])
