#pip freeze > requirements.txt

from hw2skeleton import *
from hw2skeleton import io
import sys
sys.path.append('../hw2skeleton/')
from hw2skeleton import cluster
import os


print('testing!')
similarity_metric('./data/')
