#!/usr/bin/env python

import os
import re
import sys

class BoxData:
    def __init__(self):
        self.heading_line = None
        self.level_steps = None
        self.box_bounds = []
        self.level_line = None
        self.nboxes = -1

def doit(header_file):

    # read in the Header file

    with open(header_file) as hf:
        lines = hf.readlines()

    # rename the Header to a backup

    #os.rename(header_file, header_file + ".backup")

    # read in the fields of the current header.  For data that is not
    # going to be changed, we leave it as a string so it is easy to
    # write back out

    # first line is description
    desc = lines.pop(0)

    # next comes the number of variables
    nvar = int(lines.pop(0))

    # now the variable names
    vars = []
    for _ in range(nvar):
        vars.append(lines.pop(0))

    # now the dimensionality
    ndim = int(lines.pop(0))

    # next the simulation time
    time = lines.pop(0)

    # and now the maximum level
    max_level = int(lines.pop(0))

    # next come the physical domain lo and high
    domlo = lines.pop(0)
    domhi = lines.pop(0)

    # next the set of ref_ratios -- there will be levels - 1 of these.
    # If we have max_level = 0, then this line is empty in the Header
    ref_ratios = lines.pop(0).strip().split()

    # now the domain information -- this will take the form of a group like
    # ((0,0,0) (1023,1023,127) (0,0,0)) for each level.  We'll split in via regex
    domain_string = lines.pop(0)
    domain_size = re.findall("\(\(.*?\)\)", domain_string)

    # step info (per level)
    steps = lines.pop(0).strip().split()

    # now the cell size info, on a level basis
    cell_size = []
    for _ in range(max_level+1):
        cell_size.append(lines.pop(0))

    # coordinate system
    coord = lines.pop(0)

    # finally a 0
    zero = lines.pop(0).strip()
    if not zero == "0":
        sys.exit("something went wrong")

    box_info = []

    # finally we have the box info
    for n in range(max_level+1):
        print(f"n = {n}")

        boxes = BoxData()

        # first a header that gives: level, number of boxes, and time
        tmp = lines.pop(0)
        boxes.heading_line = tmp
        level, nboxes, _ = tmp.strip().split()

        boxes.level_steps = lines.pop(0)

        # safety check
        if not int(level) == n:
            sys.exit("Error reading the box info")

        nboxes = int(nboxes)
        boxes.nboxes = nboxes

        # now we read all the box info.  For each box, there are 3 ndim lines
        for _ in range(nboxes):
            for _ in range(ndim):
                boxes.box_bounds.append(lines.pop(0))

        # now we should have a line that is of the form Level0/Cell
        boxes.level_line = lines.pop(0)
        print(boxes.level_line)

        box_info.append(boxes)

    # now output back with one fewer level
    with open("Header.tmp", "w") as hout:
        hout.write(desc)
        hout.write(f"{nvar}\n")
        for v in vars:
            hout.write(v)
        hout.write(f"{ndim}\n")
        hout.write(time)

        # we are removing a level
        hout.write(f"{max_level-1}\n")

        hout.write(domlo)
        hout.write(domhi)

        # remove one of the ref_ratios
        ref_ratio_str = " ".join([ref_ratios[k] for k in range(max_level-1)])
        hout.write(ref_ratio_str + "\n")

        # domain info -- cut off the last
        domain_size_str = " ".join([domain_size[k] for k in range(max_level)])
        hout.write(domain_size_str + "\n")

        # now the steps -- cut off the last
        steps_str = " ".join([steps[k] for k in range(max_level)])
        hout.write(steps_str + "\n")

        # now the cell dimensions -- only print out max_level-1 lines
        for n in range(max_level):
            hout.write(cell_size[n])

        # now the coord info
        hout.write(coord)

        # and the 0
        hout.write("0\n")

        # now the box data for the levels we care about
        for n in range(max_level):
            hout.write(box_info[n].heading_line)
            hout.write(box_info[n].level_steps)
            for _ in range(box_info[n].nboxes):
                for _ in range(ndim):
                    hout.write(box_info[n].box_bounds.pop(0))
            hout.write(box_info[n].level_line)


if __name__ == "__main__":

    if len(sys.argv) == 1:
        sys.exit("need to specify the plotfile on the command line")

    plot_file = sys.argv[-1]

    header = os.path.normpath(plot_file) + "/Header"

    if not os.path.isfile(header):
        sys.exit("Header not found")

    doit(header)
