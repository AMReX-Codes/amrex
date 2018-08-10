from .Bucket import Bucket

class Constants(Bucket):
    """
    The purpose of this class is to be hold user defined constants
    The constants will be concatenated into names and values string.
    """
    def __init__(self):
        Bucket.__init__(self, 'constants')

    def __setattr__(self, name, value):
        # Make sure that any constants redefined have a consistent value
        if name in self.argvattrs:
            assert self.argvattrs[name] == value, Exception('In consistent values given for user defined constants')
        Bucket.__setattr__(self, name, value)

    def attrlist(self):
        "Concatenate the attributes into a string"
        if self.argvattrs:
            names = ''
            values = ''
            for attr, value in self.argvattrs.items():
                names += ' ' + attr
                values += ' {}'.format(value)
            return ['constants.use_my_constants = 1',
                    'constants.constant_names = ' + names,
                    'constants.constant_values = ' + values]
        else:
            return []


constants = Constants()
