# Simple script to make a CAMI "profile" using the output results of CMash.
# Not a real profile since CMash doesn't give relative abundances, but will let us use OPAL to get TP/FP/etc.
import pandas as pd
import argparse
import sys
import os


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="This script puts the CMash outputs into a profile-like format that is suitable for use by OPAL for "
					"assessment purposes. IMPORTANT: IGNORE THE RELATIVE ABUNDANCES", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-t', '--taxonomy_file', type=str, help="File name of taxonomy")
	parser.add_argument('-c', '--cmash_file', type=str, help="File name of Cmash output")
	parser.add_argument('-o', '--out_file', type=str, help="File name of output *.profile file")
	parser.add_argument('--threshold', type=float, help="Threshold for which to include results in profile file", default=0)
	parser.add_argument('--key', type=str, help="Which k-mer size to do the thresholding. eg. 'k=25'. Default is largest k-mer size", default=False)
	parser.add_argument('--sample_id', type=str, help="Sample ID", default="")

	# read in the arguments
	args = parser.parse_args()
	name_and_taxpath_file = args.taxonomy_file  # '/home/dkoslicki/Data/MiCOPMinHash/Test/File_name_and_taxonomy.txt'
	cmash_out_file = args.cmash_file  # '/home/dkoslicki/Data/MiCOPMinHash/Test/Test.csv'
	profile_out_file = args.out_file  # '/home/dkoslicki/Data/MiCOPMinHash/Test/Test.profile'
	coverage_threshold = args.threshold
	sort_key = args.key
	sample_id = args.sample_id  #os.path.abspath(cmash_out_file)

	abbrv_to_rank = dict()
	abbrv_to_rank['k'] = "superkingdom"
	abbrv_to_rank['p'] = "phylum"
	abbrv_to_rank['c'] = "class"
	abbrv_to_rank['o'] = "order"
	abbrv_to_rank['f'] = "family"
	abbrv_to_rank['g'] = "genus"
	abbrv_to_rank['s'] = "species"
	abbrv_to_rank['t'] = "strain"

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

	df = pd.read_csv(cmash_out_file, index_col=0)
	max_key = df.keys()[-1]
	if not sort_key:
		sort_key = max_key
	cmash_thresh = df[df[sort_key] > coverage_threshold].sort_values(max_key, ascending=False)
	names_passed_thresh = list(cmash_thresh.index)

	fid = open(profile_out_file, 'w')
	fid.write("""#CAMI Submission for Taxonomic Profiling
@Version:0.9.1
@SampleID:%s
@Ranks:superkingdom|phylum|class|order|family|genus|species|strain

@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n""" % sample_id)

	tax_id_to_abund = dict()
	for name in names_passed_thresh:
		if name in name_to_taxpath:
			tax_info = name_to_taxpath[name]
			# add all higher ranks as well
			for rank_ind in range(1, len(tax_info[2].split('|')) + 1):
				tax_info_to_rank = tax_info[2].split('|')[:rank_ind]  # tax info up to the rank under consideration
				tax_id = [i.split('_')[2] for i in tax_info_to_rank][-1]  # get the last guy's tax id
				rank = abbrv_to_rank[tax_info_to_rank[-1].split('_')[0]]  # get the last guy's rank
				tax_path_list = [i.split('_')[2] for i in tax_info_to_rank]
				if len(tax_path_list) != len(set(tax_path_list)):
					continue
				tax_path = '|'.join([i.split('_')[2] for i in tax_info_to_rank])  # join up the tax path
				tax_path_sn = '|'.join([' '.join(i.split('_')[3:]) for i in tax_info_to_rank])  # join up the tax path names
				if tax_id in tax_id_to_abund:  # Don't care about the abundances
					tax_id_to_abund[tax_id] = .01  # so don't increment this
				else:
					tax_id_to_abund[tax_id] = .01
				percentage = tax_id_to_abund[tax_id]
			# this was the old way when I was just adding species
			#tax_id = tax_info[1]
			#rank = abbrv_to_rank[tax_info[2].split('|')[-1].split('_')[0]]
			#tax_path = '|'.join([i.split('_')[2] for i in name_to_taxpath[name][2].split('|')])
			#tax_path_sn = '|'.join([' '.join(i.split('_')[3:]) for i in name_to_taxpath[name][2].split('|')])
			#percentage = 1
				fid.write("%s\t%s\t%s\t%s\t%f\n" % (tax_id, rank, tax_path, tax_path_sn, percentage))
		else:
			print("uh oh, this file isn't in the taxonomy: %s" % name)
			sys.exit(1)
	fid.close()

# then use something like the following:
# import ProfilingTools as PF
# profile=PF.Profile('/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.profile')
# profile.normalize()
# profile._populate_missing_dont_use()
# profile.normalize()
# profile.write_file('/home/dkoslicki/Data/MiCOPMinHash/Test/RL_S001__insert_270_classified.profile.norm')
