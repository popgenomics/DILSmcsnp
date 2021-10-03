#!/usr/bin/python
import os
import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')

for line in infile:
	outfile.write(line)

infile.close()
outfile.close()

