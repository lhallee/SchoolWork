import pandas as pd
from tqdm import tqdm
from pyfaidx import Fasta
from Bio import SeqIO


def extract_id(header):  # for getting uniprot sequences
    return header.split('|')[1]


def translate(dna_seq):  # This function translates dna sequences to amino acids, where - is a stop
    dna_seq = dna_seq.upper()
    aa_seq = ''
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '-', 'TAG': '-',
        'TGC': 'C', 'TGT': 'C', 'TGA': '-', 'TGG': 'W',
    }
    if len(dna_seq) % 3 == 0:
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i + 3]
            if codon_table[codon] == '-':
                break
            aa_seq += codon_table[codon]

    else:  # not multiple of three, cannot translate
        return '0'
    return aa_seq


def mismatch(seqa, seqb):
    mismatches = []
    for i, (aminoa, aminob) in enumerate(zip(seqa, seqb)):
        if aminoa != aminob:
            mismatches.append([aminoa, aminob, i])
    return mismatches

ccds_to_uniprot = pd.read_csv('CCDS2UniProtKB.current.txt', delimiter='\t')  # load ccds ids
ccds_to_uni_dict = dict(zip(ccds_to_uniprot['#ccds'], ccds_to_uniprot['UniProtKB']))  # to dictionary
ccds_sequences = SeqIO.parse(open('CCDS_nucleotide.current.fna'),'fasta')  # load ccds sequences
uni_sequences = Fasta('uniprot_sprot.fasta', key_function=extract_id)  # load uniprot sequences

total, wrong = 0, 0
list_of_mismatch, list_of_mismatch_seqs = [], []
for fasta in tqdm(ccds_sequences):
    try:
        name, dna_seq = fasta.id, str(fasta.seq)  # get id and sequence
        ccds_id = name.split('|')[0].replace('>', '')  # filter out ccds id
        uniprot_id = ccds_to_uni_dict[ccds_id]
        prot_seq = str(uni_sequences[uniprot_id])  # get uniprot sequence of ccds id
        if prot_seq != translate(dna_seq):
            wrong += 1
            list_of_mismatch.append([ccds_id, uniprot_id])
            list_of_mismatch_seqs.append([translate(dna_seq), prot_seq])
        total += 1
    except:
        continue

print(total, wrong)


for i in range(len(list_of_mismatch_seqs)):
    ccds = list_of_mismatch_seqs[i][0]
    uni = list_of_mismatch_seqs[i][1]
    mismatches = mismatch(ccds, uni)
    print(list_of_mismatch[i][0] + '\t' + list_of_mismatch[i][1] + '\n'
            'CCDS:    ' + ccds + '\n'
            'UniProt: ' + uni + '\n',
            len(mismatches), mismatches)









