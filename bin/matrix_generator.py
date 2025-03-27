#!/usr/bin/env python

import sys
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

sample=sys.argv[1]
work_dir=sys.argv[2]
genome=sys.argv[3]

matrices = matGen.SigProfilerMatrixGeneratorFunc(sample, genome, work_dir)

