# workaround for python2/3 compatibility without six
try:
  basestring
except NameError:
  basestring = str

isstr = lambda s: isinstance(s, basestring)
