with open("particles", "r") as f:
    f.readline()
    for line in f:
        vals = line.strip().split()
        assert(vals[0] == vals[5])
        assert(vals[1] == vals[6])

print("pass!")
