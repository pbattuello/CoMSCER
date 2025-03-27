#!/usr/bin/env python3

import sys

sample=sys.argv[1]
path_script=sys.argv[2]
output_dir=sys.argv[3]
condition=sys.argv[4]
genome=sys.argv[5]

cosmic2=path_script+"/cosmic_ref/COSMIC_v2_SBS_"+genome+".txt"
cosmic3_2=path_script+"/cosmic_ref/COSMIC_v3.2_SBS_"+genome+".txt"

from SigProfilerAssignment import Analyzer as Analyze

Analyze.cosmic_fit(samples=sample,
	output=output_dir + "/TMP/" + condition + "/C2",
    genome_build=genome,
	signature_database=cosmic2)

Analyze.cosmic_fit(samples=sample,
        output=output_dir + "/TMP/" + condition + "/C3.2",
        genome_build=genome,
        signature_database=cosmic3_2)
	
