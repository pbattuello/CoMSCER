#!/usr/bin/env python

from SigProfilerMatrixGenerator import install as genInstall

print("GRCh38 installing...")
genInstall.install('GRCh38', rsync=False, bash=True)
print("done")
print("GRCh37 installing...")
genInstall.install('GRCh37', rsync=False, bash=True)
print("done")
print("Installation Completed!")
