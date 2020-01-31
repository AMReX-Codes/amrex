# Copyright 2019-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

lasers = Bucket('lasers', nlasers=0, names=[])
lasers_list = []

def newlaser(name):
    result = Bucket(name)
    lasers_list.append(result)
    lasers.nlasers += 1
    lasers.names.append(name)
    return result
