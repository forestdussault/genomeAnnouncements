#! /usr/env/python
import os
__author__ = 'adamkoziol'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    import errno
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise

# Initialise globalcount
globalcount = 0


def globalcounter():
    """Resets the globalcount to 0"""
    global globalcount
    globalcount = 0


def dotter():
    """Prints formatted time to stdout at the start of a line, as well as a "."
    whenever the length of the line is equal or lesser than 80 "." long"""
    # import time
    import sys
    # Use a global variable
    global globalcount
    if globalcount <= 80:
        sys.stdout.write('.')
        globalcount += 1
    else:
        # sys.stdout.write('\n[%s] .' % (time.strftime("%H:%M:%S")))
        sys.stdout.write('\n.')
        globalcount = 1


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = x if x else {}
        # start = (lambda y: y if y else {})(x)
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        if value:
            self.datastore[key] = value
        elif value is False:
            self.datastore[key] = value
        else:
            self.datastore[key] = "NA"


class MetadataObject(object):
    """Object to store static variables"""
    def __init__(self):
        """Create datastore attr with empty dict"""
        super(MetadataObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        """:key is retrieved from datastore if exists, for nested attr recursively :self.__setattr__"""
        if key not in self.datastore:
            self.__setattr__(key)
        return self.datastore[key]

    def __setattr__(self, key, value=GenObject(), **args):
        """Add :value to :key in datastore or create GenObject for nested attr"""
        if args:
            self.datastore[key].value = args
        else:
            self.datastore[key] = value

    def dump(self):
        """Prints only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = {}
        for attr in self.datastore:
            metadata[attr] = {}
            if not attr.startswith('__'):
                if isinstance(self.datastore[attr], str):
                    metadata[attr] = self.datastore[attr]
                else:
                    metadata[attr] = self.datastore[attr].datastore
        return metadata
