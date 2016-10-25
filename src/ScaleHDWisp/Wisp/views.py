from django.shortcuts import render
from django.http import HttpResponse

# Create your views here.
def index(request):
    context_dict = {'boldmessage': "one, two, three, four"}
    return render(request, 'Wisp/index.html', context=context_dict)
