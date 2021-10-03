#!/usr/bin/python
import sys

timeStamp = sys.argv[1]
projectpath = sys.argv[2]

models = ['SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N']
#ss = c(1, 4:11, 14:17, 19:22, 24:31, 40:41) # model_comp_2pop_all_Models_gw.R

outfile = open('{0}/{1}/{1}/gof.txt'.format(projectpath, timeStamp), 'w')

# observed data
infile = open('{0}/{1}/{1}/ABCstat_global.txt'.format(projectpath, timeStamp), 'r')

line = infile.readline().strip()
line += '\torigin\n'
outfile.write(line)

line = infile.readline().strip()
line += '\tobserved\n'
outfile.write(line)
infile.close()

# simulated data 
for model in models:
	infile = open('{0}/{1}/{1}/modelComp/{2}_0/ABCstat.txt'.format(projectpath, timeStamp, model), 'r')
	line = infile.readline().strip()
	for line in infile:
		line = line.strip()
		line += '\t{0}\n'.format(model)
		outfile.write(line)
	infile.close()
outfile.close()

