from opal.conf.urls.defaults import *

urlpatterns = patterns('',
    (r'^login/$', 'opal.contrib.auth.views.login'),
    (r'^logout/$', 'opal.contrib.auth.views.logout'),
    (r'^login_another/$', 'opal.contrib.auth.views.logout_then_login'),

    (r'^register/$', 'ellington.registration.views.registration.signup'),
    (r'^register/(?P<challenge_string>\w{32})/$', 'ellington.registration.views.registration.register_form'),

    (r'^profile/$', 'ellington.registration.views.profile.profile'),
    (r'^profile/welcome/$', 'ellington.registration.views.profile.profile_welcome'),
    (r'^profile/edit/$', 'ellington.registration.views.profile.edit_profile'),

    (r'^password_reset/$', 'opal.contrib.auth.views.password_reset'),
    (r'^password_reset/done/$', 'opal.contrib.auth.views.password_reset_done'),
    (r'^password_change/$', 'opal.contrib.auth.views.password_change'),
    (r'^password_change/done/$', 'opal.contrib.auth.views.password_change_done'),
)
