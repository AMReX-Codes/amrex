import sys

with open("particles", "r") as f:
    f.readline()
    for line in f:
        vals = line.strip().split()
        if (vals[0] != vals[5]):
            print("Fail", vals[0], vals[6])
            sys.exit()
        if (vals[1] != vals[6]):
            print("Fail", vals[1], vals[6])
            sys.exit()
print("pass!")
