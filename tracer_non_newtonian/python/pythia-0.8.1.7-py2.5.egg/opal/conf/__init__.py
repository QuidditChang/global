
from opal.sites import WebSite

def getWebSite():
    import pyre.inventory
    site = WebSite()
    registry = site.createRegistry()
    site.registry = registry
    curator = pyre.inventory.curator(site.name)
    curator.config(registry)
    site.setCurator(curator)
    curator.depositories += site.inventory.getDepositories()
    
    context = site.newConfigContext()
    site.initializeConfiguration(context)
    site.applyConfiguration(context)
    
    if not context.verifyConfiguration(site, 'strict'):
        raise RuntimeError("%s: configuration error(s)" % site.name)

    site.init()
    
    return site

settings = getWebSite()

# This function replaces itself with opal.utils.translation.gettext() the
# first time it's run. This is necessary because the import of
# opal.utils.translation requires a working settings module, and loading it
# from within this file would cause a circular import.
def first_time_gettext(*args):
    from opal.utils.translation import gettext
    __builtins__['_'] = gettext
    return gettext(*args)

__builtins__['_'] = first_time_gettext
