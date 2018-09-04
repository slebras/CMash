# simple script to get the TP's and FP's
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

name_and_taxpath_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/File_name_and_taxonomy.txt'
# low
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.csv'
#truth_file = '/home/dkoslicki/Dropbox/Repositories/firstchallenge_evaluation_Backup/profiling/MyAnalysis/GroundTruth/all/CAMI:low:pool.profile'

# medium
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RM_S001_classified.csv'
#truth_file = '/home/dkoslicki/Dropbox/Repositories/firstchallenge_evaluation_Backup/profiling/MyAnalysis/GroundTruth/all/CAMI:medium:1.profile'

# high
#cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RH_S001_classified.csv'
#truth_file = '/home/dkoslicki/Dropbox/Repositories/firstchallenge_evaluation_Backup/profiling/MyAnalysis/GroundTruth/all/CAMI:high:1.profile'

# reduced, low
cmash_out_file = '/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001_reduced_classified.txt'
truth_file = '/home/dkoslicki/Dropbox/Repositories/firstchallenge_evaluation_Backup/profiling/MyAnalysis/GroundTruth/all/CAMI:low:pool.profile'

coverage_threshold = 0
sort_key = 'k=60'  # 'k=30' barely adds any TP's with a huge hit of lots of FP's

abbrv_to_rank = dict()
abbrv_to_rank['k'] = "superkingdom"
abbrv_to_rank['p'] = "phylum"
abbrv_to_rank['c'] = "class"
abbrv_to_rank['o'] = "order"
abbrv_to_rank['f'] = "family"
abbrv_to_rank['g'] = "genus"
abbrv_to_rank['s'] = "species"
abbrv_to_rank['t'] = "strain"

rank_to_abbrv = dict()
for key, value in abbrv_to_rank.items():
	rank_to_abbrv[value] = key

# read in the taxonomy
fid = open(name_and_taxpath_file, 'r')
name_to_taxpath = dict()
for line in fid.readlines():
	line = line.strip()
	ref_file = line.split()[0]  # first column is file name
	tax_info = line.split()[1:]  # rest is the taxpath info
	if ref_file in name_to_taxpath:
		print("Uh oh, file names are not unique! culprit: %s" % ref_file)
		sys.exit(1)
	else:
		name_to_taxpath[ref_file] = tax_info
fid.close()

# Now read in the ground truth
true_taxa = dict()
for key in abbrv_to_rank.keys():
	true_taxa[key] = set()

fid = open(truth_file, 'r')
for line in fid.readlines():
	line = line.strip()
	if line and line[0] != '@':
		split_line = line.split()
		tax_id = split_line[0]
		rank = split_line[1]
		rank_abbrv = rank_to_abbrv[rank]
		true_taxa[rank_abbrv].add(tax_id)
fid.close()


TP_vers_cov = []
FP_vers_cov = []
FN_vers_cov = []
df = pd.read_csv(cmash_out_file, index_col=0)
max_key = df.keys()[-1]

# loop over coverages, get binary stats with that threshold
cov_range = np.linspace(.6, 0, 50)
#cov_range = [0.05]  # for the low complexity sample: TP=18, FP=101, FN=5
#cov_range = [0.004]  # good for the medium complexity sample: TP=54, FP=536, FN=18
#cov_range = [0.008]  # good for the high complexity sample: TP=153, FP=1123, FN=90
for coverage_threshold in cov_range:
	if not sort_key:
		sort_key = max_key
	cmash_thresh = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)
	names_passed_thresh = list(cmash_thresh.index)

	pred_taxa = dict()
	for key in abbrv_to_rank.keys():
		pred_taxa[key] = set()

	for name in names_passed_thresh:
		if name in name_to_taxpath:
			tax_info = name_to_taxpath[name]
			# add all higher ranks as well
			for rank_ind in range(1, len(tax_info[2].split('|')) + 1):
				tax_info_to_rank = tax_info[2].split('|')[:rank_ind]  # tax info up to the rank under consideration
				tax_id = [i.split('_')[2] for i in tax_info_to_rank][-1]  # get the last guy's tax id
				rank_abbrv = tax_info_to_rank[-1].split('_')[0]
				tax_path_list = [i.split('_')[2] for i in tax_info_to_rank]
				pred_taxa[rank_abbrv].update(tax_path_list)

	num_TPs = []
	num_FPs = []
	num_FNs = []
	TPs = []
	FPs = []
	FNs = []
	for rank_abbrv in 'kpcofgs':
		TP = true_taxa[rank_abbrv].intersection(pred_taxa[rank_abbrv])
		FP = pred_taxa[rank_abbrv] - true_taxa[rank_abbrv]
		FN = true_taxa[rank_abbrv] - pred_taxa[rank_abbrv]
		TPs.append(TP)
		FPs.append(FP)
		FNs.append(FN)
		num_TPs.append(len(TP))
		num_FPs.append(len(FP))
		num_FNs.append(len(FN))

	TP_vers_cov.append(num_TPs[-1])
	FP_vers_cov.append(num_FPs[-1])
	FN_vers_cov.append(num_FNs[-1])

together = np.array([TP_vers_cov, FP_vers_cov]).transpose()
plt.plot(cov_range, together)
plt.show()