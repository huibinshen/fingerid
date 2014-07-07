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

sys.path.append("../../fingerid_1.4") # path to fingerid package         
from fingerid.preprocess.msparser import MSParser
from fingerid.preprocess.fgtreeparser import FragTreeParser
from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel
from fingerid.kernel.fgtreekernel import FragTreeKernel
from fingerid.kernel.mskernel import Kernel
from fingerid.kernel.mkl import mkl
from fingerid.model.internalCV_mp import internalCV_mp
from fingerid.model.internalCV import internalCV

from fingerid.model.trainSVM import trainModels
from fingerid.model.predSVM import predModels
from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel

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
        train_km = kernel.compute_kernel(trees,km_type)
        kernel.write_kernel(train_km, km_f)
    print "Writing %s kernel to %s" % (km_type, km_f)


def kernelMP(types, ms_folder, fgtree_folder):
    """ Compute kernel with multiprocess """
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

def trainSVM(km_f, labels_f, np = 4, c_sel=False):
    """ Train svm with the specified kernel """
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
    cvpreds = internalCV_mp(train_km, labels, 5, select_c=c_sel, n_p=np)
    numpy.savetxt(cvpred_f, cvpreds, fmt="%d")
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
    types = ["PPK","NB","NI","LB","LC","LI","RLB","RLI","CPC","CP2","CPK","CSC"]
    # km_fs is the list of file names of all the kernel matrices
    km_fs = kernelMP(types, ms_folder, fgtree_folder) 
    print

    # Predict for individual kernel, using np processes
    # In the training, doing selection of best C in SVM.
    for km_f in km_fs:        
        trainSVM(km_f, fingerprints, np=4, c_sel=True)
    print

    # Get the combined kernel using 'UNIMKL', 'ALIGN' or 'ALIGNF'
    # 'ALIGNF' needs package cvxopt installed
    km_list = []
    for km_f in km_fs:
        km_list.append(numpy.loadtxt(km_f))
    output = numpy.loadtxt(fingerprints)
    ckm, w = mkl(km_list, output, 'ALIGN')
    print "Kernel weights:"
    print w
    ckm_f = 'ALIGN_kernel.txt'
    numpy.savetxt(ckm_f, ckm)
    # predict for combined kernel, using np processes
    trainSVM(ckm_f, fingerprints,  np=4, c_sel=True)


