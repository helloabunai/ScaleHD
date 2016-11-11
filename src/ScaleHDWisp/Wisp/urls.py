from django.conf.urls import url
from Wisp import views
urlpatterns = [url(r'^$', views.index, name='index'),
               ]