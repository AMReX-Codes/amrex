#!/usr/bin/env python3

import sys, os, glob, operator

if __name__ == "__main__":
    dt = float(sys.argv[3])-float(sys.argv[2])
    hours, rem = divmod(dt, 3600)
    minutes, seconds = divmod(rem, 60)
    dtstr = str(int(seconds)) + " seconds"
    if minutes > 0:
        dtstr = str(int(minutes)) + " minutes " + dtstr
    if hours > 0:
        dtstr = str(int(hours)) + " hours " + dtstr
    print("Total build time is", dtstr)
    print("More details are available at", sys.argv[1])
    log_file_name = sys.argv[1]
    log_file_dir = os.path.dirname(log_file_name)
    log_files = glob.glob(os.path.join(log_file_dir,"*.log"))
    build_time_results = {}
    for logf in log_files:
        f = open(logf,'r')
        t0 = float(f.readline())
        t1 = float(f.readline())
        build_time_results[os.path.basename(logf)[:-4]] = t1-t0
        f.close()
    f = open(log_file_name,'w')
    f.write("# (File Name, Built Time in seconds)\n")
    for it in sorted(build_time_results.items(), key=operator.itemgetter(1),reverse=True):
        f.write(str(it)+'\n')
    f.close()
