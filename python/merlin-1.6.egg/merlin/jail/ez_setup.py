
import os, sys

DEFAULT_VERSION = "0.6c3"
DEFAULT_URL     = "http://cheeseshop.python.org/packages/%s/s/setuptools/" % sys.version[:3]

def use_setuptools(version=DEFAULT_VERSION,
                   download_base=DEFAULT_URL,
                   to_dir=os.curdir,
                   download_delay=15):
    return
