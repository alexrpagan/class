import sys
import os

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def validdir(path, errmsg):
    if path is not None and os.path.isdir(path):
        return True
    else:
        sys.stderr.write(errmsg + '\n')
        return False
