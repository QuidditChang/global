def __bootstrap__():
   global __bootstrap__, __loader__, __file__
   import sys, merlin, imp
   __file__ = merlin.resource_filename(__name__,'_namemapper.so')
   del __bootstrap__, __loader__
   imp.load_dynamic(__name__,__file__)
__bootstrap__()
