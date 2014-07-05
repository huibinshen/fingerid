"""
=============================
Corss validation on the data.
=============================

"""


from svmutil import *
import commands
import numpy
import random

def internalCV(kernel, labels, n_folds, pred_f, select_c=False):
   """
    Internel cross validation using train data.

    Parameters:
    -----------
    kernel, numpy 2d array, n_train*n_train, training kernel

    labels, numpy 2d array, n_train*n_fingerprints, training labels

    n_folds, number of cross validations

    pred_f, string, file to store predicted fingerprints in the CV

    select_c, bool, whether doing C selection in CV

    Note:

    Wrtie the cross validation predition fingerprints in pred_f
    """
    (n_x, n_x) = kernel.shape
    (n_x, n_y) = labels.shape
    x = kernel
    #cv_accs = []
    pred_fp = numpy.zeros((n_x, n_y))

    tags = _label_folds(n_x, n_folds)

    for i in range(n_y):
        y = labels[:,i]
        if select_c:
            pred_fp_i, cv_acc = _CV_BestC(x, y, i, tags, n_folds)
            pred_fp[:,i] = pred_fp_i
    #        cv_accs.append(cv_acc)
        else:
            pred_fp_i, cv_acc = _CV(x, y, i, tags, n_folds)
            pred_fp[:,i] = pred_fp_i
    #        cv_accs.append(cv_acc)

    #w = open(cv_acc_f,"w")
    #w.write(" ".join(map(str,cv_accs)))
    #w.close()
    numpy.savetxt(pred_f, pred_fp, fmt="%d")

def _CV(x, y, ind, tags, n_folds):
    """
    Internel cross validation using c = 1
    """
    #print "cv on %d'th fingerprint" % ind
    n = len(x)
    pred = numpy.zeros(n)
    for i in range(1,n_folds+1):
        test = tags == i
        train = ~(tags == i)
        test = numpy.array(range(n))[test].tolist()
        train = numpy.array(range(n))[train].tolist()

        train_km = x[numpy.ix_(train,train)]
        test_km = x[numpy.ix_(test,train)]
        train_label = y[train,:]
        test_label = y[test,:]
        n_train = len(train_km)
        n_test = len(test_km)

        train_km = numpy.append(numpy.array(range(1,n_train+1)).reshape(n_train,1), train_km,1).tolist()
        test_km = numpy.append(numpy.array(range(1,n_test+1)).reshape(n_test,1), test_km,1).tolist()

        prob = svm_problem(train_label, train_km, isKernel=True)
        param = svm_parameter('-t 4 -c 1 -b 0 -q')
        m = svm_train(prob,param)
        p_label, p_acc, p_val=svm_predict(test_label,test_km, m,'-b 0 -q')
        pred[numpy.ix_(test)] = p_label

    acc = sum(pred == y) / float(n)
    return pred, acc

def _CV_BestC(x, y, ind, tags, n_folds):
    """
    Internel cross validation using c = 1
    """
    #print "cv on %d'th fingerprint" % ind
    n = len(x)

    C = numpy.array([2**-5,2**-4,2**-3,2**-2,2**-1,2**0,2**1,2**2,2**3,
     2**4,2**5,2**6,2**7,2**8,2**9,2**10])
    n_c = len(C)
    pred = numpy.zeros((n_c,n))
    for i in range(1,n_folds+1):
        test = tags == i
        train = ~(tags == i)
        test = numpy.array(range(n))[test].tolist()
        train = numpy.array(range(n))[train].tolist()

        train_km = x[numpy.ix_(train,train)]
        test_km = x[numpy.ix_(test,train)]
        train_label = y[train,:]
        test_label = y[test,:]
        n_train = len(train_km)
        n_test = len(test_km)

        train_km = numpy.append(numpy.array(range(1,n_train+1)).reshape(n_train,1), train_km,1).tolist()
        test_km = numpy.append(numpy.array(range(1,n_test+1)).reshape(n_test,1), test_km,1).tolist()

        for j in range(n_c):
            c = C[j]
            prob = svm_problem(train_label, train_km, isKernel=True)
            param = svm_parameter('-t 4 -c %s -b 0 -q' % c)
            m = svm_train(prob,param)
            p_label, p_acc, p_val=svm_predict(test_label,test_km, m,'-b 0 -q')
            pred[j,numpy.ix_(test)] = p_label

    accs = numpy.sum(pred == numpy.array(y),1) / float(n)
    max_inds = numpy.where(max(accs)==accs)
    ind = max_inds[0][0]
    return pred[ind,],accs[ind]

def _label_folds(n_x ,n):
    """
    labeling the n_x data by n folds. Sequential lableing.
    """
    tag = [0]*n_x
    for i in range(n_x):
        tag[i] = i%n + 1
    return numpy.array(tag)


#def _label_by_mol(spectra,n_cv):
#    """
#    Lableing the data by folds. Dividing the folds by kegg_id
#    """
#    mol_dict = {}
#    count = 1
#    for s in spectra:
#        if s.kegg_id not in mol_dict:
#            mol_dict[s.kegg_id] = count
#            count = count +1
    #print mol_dict
#    n_mol = len(mol_dict)

#    a = range(1,n_mol+1)
#    random.seed(1986)
#    random.shuffle(a)
#    count = 0
#    for cid,num in mol_dict.items():
#        mol_dict[cid] = a[count]
#        count = count +1

#    mol_ids = []
#    for s in spectra:
#        mol_ids.append(mol_dict[s.kegg_id])

#    tags = []
#    n_seg = n_mol/n_cv+1
#    for mol_num in mol_ids:
#        tags.append(mol_num/n_seg+1)
#    return numpy.array(tags)

