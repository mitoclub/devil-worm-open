import pandas as pd

sequence_id = []
species = []
starts_with = []
ends_with = []
period_size = []
copy_number = []
consensus_size = []
percent_matches = []
percent_indels = []
score = []
a = []
c = []
g = []
t = []
entropy = []
pattern = []
total_region = []

with open('../../Body/2Derived/nematoda_trf_output.dat', 'r') as data:
    for line in data.readlines():
        if line.startswith('Sequence:'):
            seq_string = line.split()
        elif line[0].isdigit() == True:
            tr_string = line.split()

            sequence_id.append(seq_string[1])
            if seq_string[2] == 'UNVERIFIED:':
                species.append(' '.join(seq_string[3:5]))
            else:
                species.append(' '.join(seq_string[2:4]))
            starts_with.append(tr_string[0])
            ends_with.append(tr_string[1])
            period_size.append(tr_string[2])
            copy_number.append(tr_string[3])
            consensus_size.append(tr_string[4])
            percent_matches.append(tr_string[5])
            percent_indels.append(tr_string[6])
            score.append(tr_string[7])
            a.append(tr_string[8])
            c.append(tr_string[9])
            g.append(tr_string[10])
            t.append(tr_string[11])
            entropy.append(tr_string[12])
            pattern.append(tr_string[13])
            total_region.append(tr_string[14])

trf_total = pd.DataFrame()
trf_total['SequenceID'] = sequence_id
trf_total['Species'] = species
trf_total['StartsWith'] = starts_with
trf_total['EndsWith'] = ends_with
trf_total['PeriodSize'] = period_size
trf_total['CopyNumber'] = copy_number
trf_total['ConsensusSize'] = consensus_size
trf_total['PercentMatches'] = percent_matches
trf_total['PercentIndels'] = percent_indels
trf_total['Score'] = score
trf_total['A'] = a
trf_total['C'] = c
trf_total['G'] = g
trf_total['T'] = t
trf_total['Entropy'] = entropy
trf_total['Pattern'] = pattern
trf_total['TotalRegion'] = total_region

trf_total.to_csv('../../Body/2Derived/12_TandemRepeatsOutput.csv', index=False)