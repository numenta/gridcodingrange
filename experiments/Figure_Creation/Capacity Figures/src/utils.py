import numpy as np
import pickle
from os import listdir
from os.path import isfile, join


def gather_data(mypath, key="width"):
    """
        Loads trials from a folder and 
        returns an array indexed by: m, k, t.
    """
    filenames = [join(mypath, f) for f in listdir(mypath) if isfile(join(mypath, f))]
    filenames = [f for f in filenames if f.endswith(".p")] 
    W = []
    for filename in filenames:
        with open(filename, 'rb') as f:
            data = pickle.load(f, encoding='latin1') 
            W.append(np.expand_dims(data[key], axis=-1) )

    W = np.concatenate(W, axis=-1)
    return W


def normalized_histogram(data, bins=200):
    h, b = np.histogram(data, bins=bins)
    h = h.astype(float)
    h= h/np.amax(h)
    return h, b[:-1]


