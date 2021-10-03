
# returns the global Fis for each species within the file
import sys
import random

samplingSize=2 # number of subsampled individuals
##
fastaFile = sys.argv[1] # name of the fasta file
region = sys.argv[2] # coding / non coding
outfileName = sys.argv[3]

fasta = open(fastaFile).readlines()
seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')


##
def sample_from_dict(d, sample=10):
	keys = random.sample(list(d), sample)
	values = [d[k] for k in keys]
	return dict(zip(keys, values))

def compute_Fis(locus, nMin=4, region="coding", bound_0=True):
	# locus: list of individuals within species. Each (individual) entry contains the 2 carried alleles
	# nMin: minimum number of individuals with accepted bases to estimate a Fis
	# list of accepted bases
#	nMin=4
	list_bases = ['A', 'T', 'G', 'C']

	# Fis = sum(num) / sum(denom)	
	num = 0.0
	denom = 0.0
	nSNPs = 0
	nLocus = 0
	
	# convert the "locus" dictionnary (usefull for Ho) into a single list (usefull for Hs)
	all_ind = []
	for ind_tmp in locus:
		for i in range(2):
			all_ind.append(locus[ind_tmp][i])
	
	# locus length
	L = len(all_ind[0])
	
	# number of individuals
	nInd = len(all_ind)

	# loop over all positions
	if region == "coding":
		start = 2
		end = L
		step = 3
	else:
		start = 0
		end = L
		step = 1
	for pos in range(start, end, step):
		# get the carried alleles/bases/nucleotides at the position "pos"
		pos_tmp = [] # records the alleles at the position pos
		# loop over individuals
		for ind_tmp in locus:
			# loop over alleles
			for allele_tmp in range(2):
				base = locus[ind_tmp][allele_tmp][pos]
				if base in list_bases:
					pos_tmp.append(base)
			
			# get the list of different alleles at this position; [A, A, A, T, A] -> [A, T]
			alleles = [ i for i in set(pos_tmp) ]
		
		# test for 1. sufficient number of accepted individuals (i.e, without N); 2. polymorphism
		if len(pos_tmp) >= nMin and len(alleles) > 1:
			nLocus = 1
			nSNPs += 1
			# get Hs = 1-sum of p**2; proportion of expected heterozygotes
			Hs = 1.0
			for i in alleles:
				Hs -= (pos_tmp.count(i)/(1.0*len(pos_tmp)))**2
			
			# get Ho = proportion of observed  heterozygotes
			nTot = 0.0
			nHeteroZ = 0.0
			for ind_tmp in locus:
				allele1 = locus[ind_tmp][0][pos]
				allele2 = locus[ind_tmp][1][pos]
				if allele1 in list_bases and allele2 in list_bases:
					nTot += 1
					if allele1 != allele2:
						nHeteroZ += 1
			Ho = nHeteroZ/nTot
			
			# num and denom
			if bound_0==True:
				if Hs>Ho:# following Hall, Nathan, et al. "Maximum likelihood estimation of individual inbreeding coefficients and null allele frequencies." Genetics research 94.3 (2012): 151-161.
					num += Hs-Ho
			else:
				num += Hs-Ho
			denom += Hs
	res = {}
	res['num'] = num
	res['denom'] = denom
	res['nSNPs'] = nSNPs
	res['nLocus'] = nLocus 
	return(res)


# list of loci
loci = [ i.split('|')[0] for i in seqName ]
loci = [ i for i in set(loci) ]

species = [ i.split('|')[1] for i in seqName ]
species = [ i for i in set(species) ]

nInd_per_species = {}
for species_tmp in species:
	nInd_per_species[species_tmp] = len(set([ j.split('|')[2] for j in seqName if j.split('|')[1]==species_tmp ]))

# get data
align = {} # align[species][locus][ind]
for i in range(len(seq)):
	locus_tmp = seqName[i].split('|')[0]
	species_tmp = seqName[i].split('|')[1]
	individual_tmp = seqName[i].split('|')[2]
	if species_tmp not in align:
		align[species_tmp] = {}
	if locus_tmp not in align[species_tmp]:
		align[species_tmp][locus_tmp] = {}
	if individual_tmp not in align[species_tmp][locus_tmp]:
		align[species_tmp][locus_tmp][individual_tmp] = []
	align[species_tmp][locus_tmp][individual_tmp].append(seq[i])

Fis = {}
for species_tmp in species:
	Fis[species_tmp] = {}
	Fis[species_tmp]['num'] = 0.0
	Fis[species_tmp]['denom'] = 0.0
	Fis[species_tmp]['nSNPs'] = 0
	Fis[species_tmp]['nLocus'] = 0
	Fis[species_tmp]['num_2ind'] = 0.0
	Fis[species_tmp]['denom_2ind'] = 0.0
	Fis[species_tmp]['nSNPs_2ind'] = 0
	Fis[species_tmp]['nLocus_2ind'] = 0
	
for species_tmp in species:
	for locus_tmp in align[species_tmp]:
		# all individuals
		res = compute_Fis(locus=align[species_tmp][locus_tmp], nMin=nInd_per_species[species_tmp]*2, region=region, bound_0=True)
		Fis[species_tmp]['num'] += res['num']
		Fis[species_tmp]['denom'] += res['denom']
		Fis[species_tmp]['nSNPs'] += res['nSNPs']
		Fis[species_tmp]['nLocus'] += res['nLocus']
		
		# 2 individuals
		res_2ind = compute_Fis(locus=sample_from_dict(d=align[species_tmp][locus_tmp], sample=samplingSize), nMin=samplingSize*2, region=region, bound_0=False)
		Fis[species_tmp]['num_2ind'] += res_2ind['num']
		Fis[species_tmp]['denom_2ind'] += res_2ind['denom']
		Fis[species_tmp]['nSNPs_2ind'] += res_2ind['nSNPs']
		Fis[species_tmp]['nLocus_2ind'] += res_2ind['nLocus']
	if Fis[species_tmp]['denom'] != 0:
		Fis[species_tmp]['Fis'] = Fis[species_tmp]['num']/Fis[species_tmp]['denom']
	else:
		Fis[species_tmp]['Fis'] = -9

	if Fis[species_tmp]['denom_2ind'] != 0:
		Fis[species_tmp]['Fis_2ind'] = Fis[species_tmp]['num_2ind']/Fis[species_tmp]['denom_2ind']
	else:
		Fis[species_tmp]['Fis_2ind'] = -9

# return results
output = "species\tnInd\tnLocus\tnSNPs\tglobal_Fis\tavg_Hs\tavg_Ho\tnLocus_2ind\tnSNPs_2ind\tglobal_Fis_2ind\n"
for i in species:
	if Fis[i]['nSNPs'] != 0: 
		output += '{species}\t{nInd}\t{nLocus}\t{nSNPs}\t{Fis:.5f}\t{Hs:.5f}\t{Ho:.5f}\t{nLocus_2ind}\t{nSNPs_2ind}\t{Fis_2ind:.5f}\n'.format(species=i, nInd=nInd_per_species[i], nLocus=Fis[i]['nLocus'], nSNPs=Fis[i]['nSNPs'], Fis=Fis[i]['Fis'], Hs=Fis[i]['denom']/Fis[i]['nSNPs'], Ho=-1*Fis[i]['num']/Fis[i]['nSNPs'] + Fis[i]['denom']/Fis[i]['nSNPs'], nLocus_2ind=Fis[i]['nLocus_2ind'], nSNPs_2ind=Fis[i]['nSNPs_2ind'], Fis_2ind=Fis[i]['Fis_2ind'])
	else:
		output += '{species}\t{nInd}\t{nLocus}\t{nSNPs}\t{Fis:.5f}\t{Hs:.5f}\t{Ho:.5f}\t{nLocus_2ind}\t{nSNPs_2ind}\t{Fis_2ind:.5f}\n'.format(species=i, nInd=nInd_per_species[i], nLocus=Fis[i]['nLocus'], nSNPs=Fis[i]['nSNPs'], Fis=Fis[i]['Fis'], Hs=-9, Ho=-9, nLocus_2ind=Fis[i]['nLocus_2ind'], nSNPs_2ind=Fis[i]['nSNPs_2ind'], Fis_2ind=Fis[i]['Fis_2ind'])


outfile = open(outfileName, 'w')
outfile.write(output)
outfile.close()


