import sys

with open("particles", "r") as f:
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    for line in f:
        vals = line.strip().split()
        if (vals[3] != vals[5] or vals[3] != vals[7]):
            print("Fail", vals[3], vals[5], vals[7])
            sys.exit()
        if (vals[4] != vals[6] or vals[4] != vals[8]):
            print("Fail", vals[4], vals[6], vals[8])
            sys.exit()
print("pass!")
