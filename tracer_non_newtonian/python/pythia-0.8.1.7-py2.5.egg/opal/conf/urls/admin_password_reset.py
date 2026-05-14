from opal.conf.urls.defaults import *

urlpatterns = patterns('opal.views',
    (r'^$', 'registration.passwords.password_reset', {'is_admin_site' : True}),
    (r'^done/$', 'registration.passwords.password_reset_done'),
)
