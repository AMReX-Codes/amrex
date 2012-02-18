#!/usr/bin/env python

import os

basename = "wd_576_base_plt"
diagname = "fthermo.Linux.PathScale.exe"

plotfiles = []

# find all the plotfiles beginning with the basename
for file in os.listdir(os.getcwd()):
    if (os.path.isdir(file) and file.startswith(basename)):
        plotfiles.append(file)


plotfiles.sort()

print plotfiles

n = 0
while (n < len(plotfiles)-1):
    print "%s --infile1 %s --infile2 %s >> cp_diag.out" % (diagname, plotfiles[n], plotfiles[n+1])
    n += 1


