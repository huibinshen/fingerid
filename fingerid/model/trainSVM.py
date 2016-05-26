"""
===================================================
Train the model using SVM on all the training data.
===================================================
"""
import sys
import commands
import pickle
import numpy
import multiprocessing

from svmutil import *

def trainModels(kernel, labels, model_dir, select_c=False, n_p=4, prob=False):
    """Train SVM on all the training data.
    All the fingeprints will be trained by n_p processes.
    Use 5 fold cross validation to find the best C to train SVM.

    Parameters:
    -----------
    kernel: numpy 2-d array, kerenl matrix

    labels: numpy 2-d array, labels of the train data; binary value, 1,-1

    model_dir, string, folder to put trained models

    select_c: bool, doing C selection or not. Default setting c = 1.

    n_p: int, number of processes to use.

    prod: boolean, set to True if want probability output.

    Note:
    -----
    The trained models will be stored in the folder MODELS

    """
    commands.getoutput("mkdir %s" % model_dir)
    print "Create directory %s to store the trained models" % model_dir
    (n_x, n_x) = kernel.shape
    (n_x, n_y) = labels.shape
    x = kernel
    # internally use 5 folds cross validation to find the best parameter
    n_folds = 5 

    tags = _label_folds(n_x, n_folds)

    #cv_accs = numpy.zeros(n_y) # cross validation accuracy
    #result_queue = multiprocessing.Queue(n_y)
    if n_y < n_p:
        print "Only %d fingerprints are used" % n_y
        print "Change n_p to %d" % n_y
        n_p = n_y
    task_dict = {}
    for i in range(n_y):
        task_dict[i%n_p] = []
    for i in range(n_y):
        task_dict[i%n_p].append(i)

    jobs = []
    for i in range(n_p):
        #y = labels[:,i]
        if select_c:
            p = multiprocessing.Process(target=_trainSVMBestC, 
                                        args=(x, labels, model_dir, 
                                              task_dict[i], tags, prob,))
            jobs.append(p)
            p.start()
        else:
            x = numpy.append(numpy.array(
                    range(1,n_x+1)).reshape(n_x,1),x,1).tolist()
            p = multiprocessing.Process(target=_trainSVM, 
                                        args=(x, labels, model_dir, 
                                              task_dict[i], prob,))
            jobs.append(p)
            p.start()

    for job in jobs:
        job.join()

    # collect result
    #for i in range(n_y):
    #    res = result_queue.get()
    #    fp_ind = res.fp_ind
        #models[fp_ind] = res.model
    #pickle.dump(models, open(model_f,"wb"))
    #return models
    #    cv_accs[fp_ind] = res.acc

    #w =open(cv_acc_f,"w")
    #w.write(" ".join(map(str,cv_accs)))
    #w.close()


#    for i in range(n_y):
#        y = labels[:,i]
#        if select_c:
#            if merge:
#                tags = _label_folds(n_x, 5)
#            else:
#                tags = _label_by_mol(spectra, 5)
#            cv_acc = _trainSVMBestC(x, y, model_dir, i, tags)
#            cv_accs.append(cv_acc)
#        else:
            # formating for libsvm
#            x = numpy.append(numpy.array(
#                    range(1,n_x+1)).reshape(n_x,1),x,1).tolist()
#            cv_acc = _trainSVM(x, y, model_dir, i)
#            cv_accs.append(cv_acc)


def _trainSVM(kernel, labels, model_dir, inds, pb):
    """
    Train the svm with c = 1.
    """
    for ind in inds:
        label = labels[:,ind]
        #rint "train %d th fingerprint now ..." % (i+1)
        prob = svm_problem(label, kernel, isKernel=True)
        if pb:
            param = svm_parameter('-t 4 -c 1 -b 1 -q')
            m = svm_train(prob, param)
            svm_save_model('%s/%d.model' % (model_dir,ind), m)
        else:
            param = svm_parameter('-t 4 -c 1 -b 0 -q')
            m = svm_train(prob, param)
            svm_save_model('%s/%d.model' % (model_dir,ind), m)
        #res = result(ind)
        #res.ind = ind
        #res.model = m
        #Queue.put(res)

def _trainSVMBestC(k_m, labels, model_dir, inds, tags, pb):
    """
    Train the svm with the best C. C is selected from 5 folds cv.
    """
    for ind in inds:
        label = labels[:,ind]
        #print "train %dth fingerprint now (select C)..." % (ind+1)
        C = numpy.array([2**-5,2**-4,2**-3,2**-2,2**-1,2**0,2**1,2**2,2**3,
             2**4,2**5,2**6,2**7,2**8,2**9,2**10])
        n = len(k_m)
        accs = []
        #res = result(ind)
        for c in C:
            pred = numpy.zeros(n) # store the predict output
            for i in range(1,6):
                test = tags==i
                train = ~(tags==i)
                test = numpy.array(range(len(k_m)))[test].tolist()
                train = numpy.array(range(len(k_m)))[train].tolist()

                train_km = k_m[numpy.ix_(train,train)]
                test_km = k_m[numpy.ix_(test,train)]
                train_label = label[train]
                test_label = label[test]
                n_train = len(train_km)
                n_test = len(test_km)

                # formating for libsvm
                train_km = numpy.append(numpy.array(range(1,n_train+1)).reshape(n_train,1), train_km,1).tolist()
                test_km = numpy.append(numpy.array(range(1,n_test+1)).reshape(n_test,1), test_km,1).tolist()

                prob = svm_problem(train_label, train_km, isKernel=True)
                param = svm_parameter('-t 4 -c %f -b 0 -q' % c)
                m = svm_train(prob,param)
                p_label, p_acc, p_val=svm_predict(test_label,test_km, m,'-b 0 -q')
                pred[numpy.ix_(test)] = p_label
            acc = sum(pred==label) / float(n)
            accs.append(acc)
        accs = numpy.array(accs)
        best_c = C[accs==max(accs)][0]

        # find the best c, using the best c to train
        kernel = numpy.append(numpy.array(
                range(1,n+1)).reshape(n,1),k_m,1).tolist()

        if pb:
            prob = svm_problem(label, kernel, isKernel=True)
            param = svm_parameter('-t 4 -c %f -b 1 -q' % best_c)
            m = svm_train(prob, param)
            svm_save_model('%s/%d.model' % (model_dir,ind), m)
        else:
            prob = svm_problem(label, kernel, isKernel=True)
            param = svm_parameter('-t 4 -c %f -b 0 -q' % best_c)
            m = svm_train(prob, param)
            svm_save_model('%s/%d.model' % (model_dir,ind), m)

        #res.ind = ind
        #res.acc = max(accs)
        #res.model = m
        #Queue.put(res)





def _label_folds(n_x ,n):
    """
    labeling the data by folds. Sequential lableing.
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
