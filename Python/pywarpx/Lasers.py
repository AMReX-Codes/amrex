from .Bucket import Bucket

lasers = Bucket('lasers', nlasers=0, names=None)
lasers_list = []

def newlaser(name):
    result = Bucket(name)
    lasers_list.append(result)
    lasers.nlasers += 1
    if lasers.names is None:
        lasers.names = name
    else:
        lasers.names += ' ' + name
    return result
