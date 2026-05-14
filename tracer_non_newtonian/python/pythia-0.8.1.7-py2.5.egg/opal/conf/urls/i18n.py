from opal.conf.urls.defaults import *

urlpatterns = patterns('',
    (r'^setlang/$', 'opal.views.i18n.set_language'),
)
