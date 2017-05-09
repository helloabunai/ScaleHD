import re
import difflib

def similar(seq1, seq2):
	return difflib.SequenceMatcher(a=seq1.lower(), b=seq2.lower()).ratio()

def scraper(intv_dict, intervening_str):
	for dna_module in re.finditer(intv_dict['Mask'], intervening_str):
		if intv_dict['Count'] == 0:
			intv_dict['StartIDX'] = dna_module.start(); intv_dict['EndIDX'] = dna_module.end()
		elif intv_dict['Count'] != 0:
			intv_dict['EndIDX'] = dna_module.end()
		intv_dict['Count'] += 1
	return intv_dict

for intervening in ['CAACAGCCGCCA','CAACAGCAACAGCCGCCA','CAACAGCCGCCACCGCCA','CAACAG','CCGCCA','CAGCAGCAGCCGCCGCCG',
					'CAACAACAGCCGCCA', 'CCTCAACAGCCGCCA', 'CAACAGCTTCCGCCA', 'CAACAGCCGCCACCTCCT', 'CATCCGCCGCCA',
					'CAACAGCTGCCA','CAACATCAGCCGCCA', 'CAACAGCCTCCACCA']:

	##
	## Set up data structure
	int_one = {'Mask': 'CAACAG', 'Count': 0, 'StartIDX': 0, 'EndIDX': 0, 'Label': '', 'Suffix': ''}
	int_two = {'Mask': 'CCGCCA', 'Count': 0, 'StartIDX': 0, 'EndIDX': 0, 'Label': '', 'Suffix': ''}
	print '\n\n>>Intervening: ', intervening

	##
	## Scrape data into structure
	for mask_dict in [int_one, int_two]:
		scraper(mask_dict, intervening)

	intervening_flag = True; atypical_flag = True; int_one_offset_flag = False; int_two_offset_flag = False
	int_one_investigate = False; int_two_investigate = False
	caacag_count = int_one['Count']; ccgcca_count = int_two['Count']
	if caacag_count == 0 and ccgcca_count == 0: intervening_flag = False
	if caacag_count != 1 and ccgcca_count != 1: atypical_flag = True
	int_one_offset = 0; int_one_simscore = 0; int_two_simscore = 0

	##########################
	##CAACAG (int one) check##
	##########################
	if int_one['Count'] > 0:
		if not int_one['StartIDX'] == 0:
			offset_str = intervening.split(int_one['Mask'])[0]
			int_one_offset = len(offset_str)
			int_one_offset_flag = True
			int_one_investigate = True
	else:
		remainder = len(intervening)%6
		if not remainder == 0:
			offset_mutated = intervening.split(intervening[remainder:remainder+6])[0]
			int_one_offset = len(offset_mutated)
			potential_mask = intervening[remainder:remainder+6]
			int_one_offset_simscore = similar('CAACAG', potential_mask)
			int_one_offset_flag = True
			if int_one_offset_simscore >= 0.5:
				int_one['Mask'] = potential_mask
				scraper(int_one, intervening)
				int_one_investigate = True
		else:
			if len(intervening) > 6 and not int_one_offset_flag:
				int_one_simscore = similar('CAACAG', intervening[0:6])
				if int_one_simscore >= 0.5:
					int_one['Mask'] = intervening[0:6]
					scraper(int_one, intervening)
					int_one_investigate = True
	print '{}: {}{}, {}{}, {}{}'.format(int_one['Mask'], 'Count ', int_one['Count'], 'Start ',int_one['StartIDX'], 'End ', int_one['EndIDX'])

	##########################
	##CCGCCA (int two) check##
	##########################
	if int_two['Count'] > 0:
		offset = (int_one['Count'] * 6) + int_one_offset
		if not int_two['StartIDX'] == offset:
			int_two_offset_flag = True
			offset_str = intervening.split(int_two['Mask'])[0].split(int_one['Mask'])[1]
			int_two_investigate = True
	else:
		remainder = len(intervening)%6
		if not remainder == 0:
			lhinge = remainder+6; rhinge = lhinge+6
			offset_mutated = intervening.split(intervening[lhinge:rhinge])[0].split(int_one['Mask'])[1]
			offset = (int_one['Count']*6)+len(offset_mutated)
			if not int_two['StartIDX'] == offset:
				offset_str = intervening[lhinge:rhinge]
				int_two_offset_simscore = similar('CCGCCA', offset_str)
				int_two_offset_flag = True
				if int_two_offset_simscore >= 0.5:
					int_two['Mask'] = offset_str
					scraper(int_two, intervening)
					int_two_investigate = True
		if len(intervening) > 6 and not int_two_offset_flag:
			int_two_simscore = similar('CCGCCA', intervening[6:12])
			if int_two_simscore >= 0.5:
				int_two['Mask'] = intervening[6:12]
				scraper(int_two, intervening)
				int_two_investigate = True
	print '{}: {}{}, {}{}, {}{}'.format(int_two['Mask'], 'Count ', int_two['Count'], 'Start ',int_two['StartIDX'], 'End ', int_two['EndIDX'])

	#################################
	##Anything longer than typical?##
	#################################
	if intervening_flag and len(intervening) > 12:
		returned_suffix = intervening[int_two['EndIDX']:]
		if returned_suffix:
			int_two_investigate = True

	###########################
	##If not present at all..##
	###########################
	if not intervening_flag:
		int_one['Count'] = 0; int_one_investigate = True
		int_two_investigate = True; int_two['Count'] = 0

	########################
	##Build genotype label##
	########################

	if int_one_investigate: int_one['Suffix'] = '*'
	if int_two_investigate: int_two['Suffix'] = '*'

	int_one['Label'] = '{}'.format(str(int_one['Count'])+int_one['Suffix'])
	int_two['Label'] = '{}'.format(str(int_two['Count'])+int_two['Suffix'])
	print '--> Genotype: {}_{}_{}_{}_{}'.format('CAG',int_one['Label'],int_two['Label'],'CCG','CCT')
