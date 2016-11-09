
class Bucket(object):
    """
    The purpose of this class is to be a named bucket for holding attributes.
    This attributes will be concatenated into a string and passed into argv during initialization.
    """
    def __init__(self, instancename):
        self._localsetattr('instancename', instancename)
        self._localsetattr('argvattrs', [])

    def _localsetattr(self, name, value):
        object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        self.argvattrs.append(name)
        #self.__dict__[name] = value
        object.__setattr__(self, name, value)

    def attrlist(self):
        "Concatenate the attributes into a string"
        result = []
        for attr in self.argvattrs:
            attrstring = '{0}.{1}={2} '.format(self.instancename, attr, getattr(self, attr))
            result += [attrstring]
        return result

