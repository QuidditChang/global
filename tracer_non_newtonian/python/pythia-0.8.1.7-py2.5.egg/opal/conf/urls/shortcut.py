from opal.conf.urls.defaults import *

urlpatterns = patterns('opal.views',
    (r'^(?P<content_type_id>\d+)/(?P<object_id>\d+)/$', 'defaults.shortcut'),
)
