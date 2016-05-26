"""
================================================================================
The fingerprint prediction on the METLIN dataset in the paper: 

Huibin Shen, Kai D\"ahrkop, Sebastian B\"ocker and Juho Rousu: Metabolite Identification through Multiple Kernel Learning on Fragmentation Trees. In the proceedings of ISMB 2014, Bioinformatics 30(12), i157-i164 (2014). 

Note that the PPK kernel improved due to tuning the variance parameters of mass 
to charge ratio (sm) and intensities (si) after the paper is published. However,
ALIGNF still reveals the best performance. 
================================================================================
"""

import sys
import numpy
import multiprocessing
import warnings; warnings.filterwarnings('ignore')
# Comment the following line if fingerid has been installed,
# otherwise, leave it there
sys.path.append("../../fingerid") # path to fingerid package         

from fingerid.preprocess.msparser import MSParser
from fingerid.preprocess.fgtreeparser import FragTreeParser
from fingerid.kernel.mskernel import Kernel
from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel
from fingerid.kernel.fgtreekernel import FragTreeKernel
from fingerid.kernel.mkl import mkl
from fingerid.model.internalCV_mp import internalCV_mp
from fingerid.model.internalCV import internalCV

def compute_kernel(km_f, ms_folder, fgtree_folder):
    """ compute the kernel indicated by keyword """

    km_type = km_f[:km_f.find("_kernel.txt")]
    print "computing kernel", km_type

    # parse ms and trees
    trees = []
    fgtreeparser = FragTreeParser()
    train_ms_list = []
    msparser = MSParser()

    train_ms_list = msparser.parse_dir(ms_folder)
    trees = fgtreeparser.parse_dir(fgtree_folder)

    if km_type == 'PPK':
        sm = 0.00001
        si = 100000
        # shoud select sm and si by cross validation
        kernel = TwoDGaussianKernel(sm, si)
        train_km = kernel.compute_train_kernel(train_ms_list)
        kernel.write_kernel(train_km, km_f)
    else:
        kernel = FragTreeKernel()
        train_km = kernel.compute_train_kernel(trees,km_type)
        kernel.write_kernel(train_km, km_f)
    print "Writing %s kernel to %s" % (km_type, km_f)

def kernelMP(types, ms_folder, fgtree_folder):
    """ Compute kernel with multiprocess. Using as many processes 
        as many trypes.
    """
    km_fs = []
    procs = []
    for km_type in types:
        km_f = km_type+ "_kernel.txt"
        km_fs.append(km_f)
        p = multiprocessing.Process(target=compute_kernel,args=(km_f,
                                    ms_folder, fgtree_folder))
        p.start()
        procs.append(p)
    for p in procs:
        p.join()    # wait until all sub-processes finished
    return km_fs

def trainSVMCV(km_f, labels_f, n_p=4, c_sel=False):
    """ Train svm with the specified kernel 
        Write the cross validations prediction of fingerprints.
    """
    km_type = km_f[:km_f.find("_kernel.txt")]
    print "Train SVM for kernel %s" % km_type
    # this files will be generated
    cvpred_f = "cvpred_" + km_type + ".txt" # the prediction file

    kernel = Kernel()
    train_km = kernel.read_kernel(km_f)
    labels = numpy.loadtxt(labels_f)

    # Use n_p process to do the cross validation on one computer.
    # For larger task, consider using computing clusters to parallel all
    # the processes.
    #cvpreds = internalCV(train_km, labels, 5, select_c=c_sel)

    prob = False # set prob = True if want probability output
    cvpreds = internalCV_mp(train_km, labels, 5, select_c=c_sel, n_p=n_p, 
                            prob=prob)
    numpy.savetxt(cvpred_f, cvpreds, fmt="%.4f")
    print "Writting prediction in %s" % cvpred_f

if __name__ == "__main__":

    # Set data info, default run on a small 50 compounds dataset.

    # To run the 978 componds dataset, change fgtree_folder to "metlin_trees".
    # The MS/MS used in the paper can be downloaded from METLIN database
    # with the same metlin id used in the filenames of fgtrees.
    fgtree_folder = "test_data/train_trees/"
    ms_folder = "test_data/train_ms/"
    fingerprints = "test_data/train_output.txt" # output we want to predict

    # compute 12 kernels using 12 process.
    # Could compute one by one by use compute_kernel function
    # Now we use as many cores as number of kernels
    # If no fragmentation trees available, just use 'PPK'
    types = ["PPK","NB","NI","LB","LC","LI","RLB","RLI","CPC","CP2","CPK","CSC"]
    # km_fs is the list of file names of all the kernel matrices
    km_fs = kernelMP(types, ms_folder, fgtree_folder) 
    print

    # Predict for individual kernel, using np processes
    # In the training, doing selection of best C in SVM.
    for km_f in km_fs:        
        trainSVMCV(km_f, fingerprints, n_p=4, c_sel=True)
    print

    # Get the combined kernel using 'UNIMKL', 'ALIGN' or 'ALIGNF'
    # 'ALIGNF' needs package cvxopt installed
    km_list = []
    for km_f in km_fs:
        km_list.append(numpy.loadtxt(km_f))
    output = numpy.loadtxt(fingerprints)

    mkl_algo = 'ALIGN' # 'ALIGNF' is the other option
    ckm, w = mkl(km_list, output, mkl_algo)
    print "Kernel weights:"
    print w
    ckm_f = '%s_kernel.txt' % mkl_algo
    numpy.savetxt(ckm_f, ckm)

    # Cross validation uses n_p processes
    # trainSVMCV function will write cross validation predictions out.
    # If you want to store the models for later use, please refer
    # train_test.py for more detail.
    trainSVMCV(ckm_f, fingerprints,  n_p=4, c_sel=True)


