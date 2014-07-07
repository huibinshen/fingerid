"""
================================================================================
Pipeline when train and test are sperated.
If only interested in cross validation on train set, use the pipeline of 
shen_ISMB2014.py.
================================================================================
"""

import sys
import numpy
import multiprocessing
import warnings; warnings.filterwarnings('ignore')

sys.path.append("../../fingerid_1.4/") # path to fingerid package         
from fingerid.preprocess.msparser import MSParser
from fingerid.preprocess.fgtreeparser import FragTreeParser
from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel
from fingerid.kernel.fgtreekernel import FragTreeKernel
from fingerid.kernel.mskernel import Kernel
from fingerid.kernel.mkl import mkl
from fingerid.model.internalCV_mp import internalCV_mp

from fingerid.model.trainSVM import trainModels
from fingerid.model.predSVM import predModels
from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel


if __name__ == "__main__":
    """ Another pipeline when you have train/test instead of cross validation""" 
    # parse data
    print "parse data\n"
    fgtreeparser = FragTreeParser()
    msparser = MSParser()
    
    train_ms = msparser.parse_dir("test_data/train_ms/")
    test_ms = msparser.parse_dir("test_data/test_ms/")
    train_trees = fgtreeparser.parse_dir("test_data/train_trees")
    test_trees = fgtreeparser.parse_dir("test_data/test_trees")
    labels = numpy.loadtxt("test_data/train_output.txt")
    n_train, n_fp = labels.shape
    n_test = len(test_trees)
    n = n_test + n_train
    all_ms = train_ms + test_ms
    all_trees = train_trees + test_trees

    # compute train and test kernels
    types = ["PPK","NB","NI","LB","LC","LI","RLB","RLI","CPC","CP2","CPK","CSC"]
    train_km_list = []
    test_km_list = []
    # can use mulitp process
    for ty in types:
        print "computing %s kernels" % ty
        if ty == "PPK":
            sm = 0.00001
            si = 100000
            # shoud select sm and si by cross validation
            kernel = TwoDGaussianKernel(sm, si)
            km = kernel.compute_train_kernel(all_ms)
            train_km = km[0:n_train, 0:n_train]
            test_km = km[n_train:n, 0:n_train]          
            train_km_list.append(train_km)
            test_km_list.append(test_km)
        else:
            kernel = FragTreeKernel()
            km = kernel.compute_kernel(all_trees,ty) 
            train_km = km[0:n_train, 0:n_train]
            test_km = km[n_train:n, 0:n_train]          
            train_km_list.append(train_km)
            test_km_list.append(test_km)

    print "combine kernels\n"
    # compute combined train and test kernel
    train_ckm, w = mkl(train_km_list, labels, 'ALIGN')
    test_ckm = numpy.zeros((n_test, n_train))
    for i in range(len(types)):
        test_ckm + w[i]*test_km_list[i]

    # train with 4 processe
    print "train models and make prediction"
    # MODELS is the folder to store trained models
    trainModels(train_ckm, labels, "MODELS", select_c=False, n_p=4)
    #print models
    preds = predModels(test_ckm, n_fp, "MODELS")
    numpy.savetxt("predictions.txt", preds, fmt="%d")








