import glob
import numpy as np
import pandas as pd
from queryuni import queryuni
from checkpdb import checkpdb
from pdb2cd import pdb2cd
from getpdbs import getpdbs
import urllib.parse
import urllib.request
from io import StringIO, BytesIO
import os

def nrmsd(spec, out, pdblist):

    #DENOMINATOR WORK
    m,n = spec.shape

    #NUMERATOR WORK
    deltamax, numerator = [], []
    for i in range(1,n):
        difference = 0
        for j in range(0,m):
            difference += (spec[j,i] - out[j,i]) ** 2
        ind_lst = list(np.fmax(spec[:,i],out[:,i]))
        ind_lst = ind_lst.index(max(ind_lst))
        deltamax.append((spec[ind_lst,i] - out[ind_lst,i]) ** 2)
        numerator.append(difference ** 0.5)

    #RMSD CALCULATION
    n = np.asfarray(numerator)
    d = np.asfarray(deltamax)

    return dict(zip(pdblist, list(n/d)))
