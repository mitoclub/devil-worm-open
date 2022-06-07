import numpy as np
import pandas as pd
from Bio import Seq
from Bio.Data import CodonTable
from Bio import SeqIO
from itertools import product
import re
import tqdm

dna_codons =[''.join(c) for c in product('ACGT', repeat = 3)]


def codon_usage(sequence: Seq, normalize = True):

    table = CodonTable.unambiguous_dna_by_id[5]
    codon_list = list(table.forward_table.keys()) + table.stop_codons
    sequence = str(sequence)
    if len(set(sequence).difference(set(['A','C', 'G', 'T'])))>0:
        return None
    start_codon_positions = [re.search(startcodon, sequence) for startcodon in table.start_codons]
    start_codon_positions = list(filter(lambda x: x is not None, start_codon_positions))
    if len(start_codon_positions) == 0:
        raise ValueError('No start codons in seq')

    startpos = np.min(list(map(lambda x:x.start(),  start_codon_positions)))

    codons = []
    translated = []
    for p in range(startpos, len(sequence), 3):
        codon = sequence[p:p+3]
        if codon not in dna_codons:
            continue

        codons.append(codon)

        if codon in table.stop_codons:
            break
        translated.append(table.forward_table[codon])

    codons, counts = np.unique(codons, return_counts=True, axis=0)
    codons = [''.join(c) for c in codons]

    ret = pd.Series(index=codon_list, data=0)
    ret[codons] = counts
    if normalize:
        ret = ret / ret.sum()
    return ret, ''.join(translated)

def cds_codon_usage(seq:Seq):
    cds_stats = []
    cds_col_names = []

    aa = []
    aa_table = pd.DataFrame()

    for feature in seq.features:
        if feature.type != 'CDS' or feature.strand != 1:
            continue
        nucl = feature.extract(seq.seq)
        cu = codon_usage(nucl, normalize=False)
        if cu is not None:
            cu, translated = cu
            cds_stats.append(cu)

            k = 'protein_id' if 'protein_id' in feature.qualifiers else 'label'
            if 'gene' in feature.qualifiers:
                cds_col_names.append(feature.qualifiers[k][0] + '|' + feature.qualifiers['gene'][0])
                aa.append(feature.qualifiers['gene'][0])
                aa.append(feature.qualifiers[k][0])
            else:
                cds_col_names.append(feature.qualifiers[k][0])
                aa.append('None')
                aa.append(feature.qualifiers[k][0])

            if translated is not None:
                aa.append(translated)
            else:
                aa.append('None')
            if 'translation' in feature.qualifiers:
                aa.append(feature.qualifiers['translation'][0])
            else:
                aa.append('None')
            aa_table = aa_table.append(pd.Series(aa), ignore_index=True)
            aa = []

    if len(cds_stats) == 0:
        return None, aa_table
    cds_stat_table = pd.concat(cds_stats, axis=1)
    cds_stat_table.columns = cds_col_names

    aa_table.columns = ['GenName', 'GenID', 'TranslatedAminoAcids', 'GenBankAminoAcids']

    return cds_stat_table, aa_table


if __name__ == '__main__':
    seq_reader = SeqIO.parse('../../Body/1Raw/mitochondrion.genomic.gbff', format='genbank')

    for sequence in tqdm.tqdm(seq_reader):
        stat_tbl, aa_tbl = cds_codon_usage(sequence)

        orgname = sequence.annotations['organism']
        if stat_tbl is not None:
            stat_tbl.to_csv(f'../../Body/2Derived/codon_usage/{sequence.name}_{orgname}.tsv', sep='\t')
        if aa_tbl is not None:
            aa_tbl.to_csv(f'../../Body/2Derived/aminoacids/{sequence.name}_{orgname}.tsv', sep='\t')
