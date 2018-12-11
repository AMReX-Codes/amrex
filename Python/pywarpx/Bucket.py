class Bucket(object):
    """
    The purpose of this class is to be a named bucket for holding attributes.
    This attributes will be concatenated into a string and passed into argv during initialization.
    """
    def __init__(self, instancename, **defaults):
        self._localsetattr('instancename', instancename)
        self._localsetattr('argvattrs', {})
        self.argvattrs.update(defaults)

    def _localsetattr(self, name, value):
        object.__setattr__(self, name, value)

    def __setattr__(self, name, value):
        self.argvattrs[name] = value

    def __getattr__(self, name):
        try:
            return self.argvattrs[name]
        except KeyError:
            return object.__getattr__(self, name)

    def check_consistency(self, vname, value, errmsg):
        if vname in self.argvattrs:
            assert (self.argvattrs[vname] is None) or (self.argvattrs[vname] == value), Exception(errmsg)

    def attrlist(self):
        "Concatenate the attributes into a string"
        result = []
        for attr, value in self.argvattrs.items():
            if value is None:
                continue
            # --- repr is applied to value so that for floats, all of the digits are included.
            # --- The strip is then needed when value is a string.
            if isinstance(value, str):
                rhs = value
            elif hasattr(value, '__iter__'):
                # --- For lists, tuples, and arrays make a space delimited string of the values
                rhs = ' '.join(map(repr, value))
            else:
                rhs = value
            attrstring = '{0}.{1} = {2}'.format(self.instancename, attr, repr(rhs).strip("'\""))
            result += [attrstring]
        return result

