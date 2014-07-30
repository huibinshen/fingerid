"""
=================================================
Predict the fingerprints using trained SVM models
=================================================
"""
from svmutil import *
import numpy

def predModels(test_km, n_fp, model_dir, prob=False):
    """
    Using trained models in model_dir and testing kernels to predict fingerprits
    
    Parameters
    ----------
    test_km: numpy.2d array, n_test * n_train. Testing kernels

    n_fp: int, number of fingerprints

    model_dir: folder of trained models

    Returns:
    -------
    pred_fp: numpy.2d array, n_test * n_fp. Predicted fingerprints matrix
             or probabilites of positive label

    """
    n_test = len(test_km)

    # formating for libsvm
    test_km = numpy.append(numpy.array(range(1,n_test+1)).reshape(n_test,1), test_km,1).tolist()

    pred_fp = numpy.zeros((n_test, n_fp))
    dummy_y = numpy.zeros(n_test)
    for i in range(n_fp):
        model_f = "%s/%d.model" % (model_dir, i)
        m = svm_load_model(model_f)
        if prob:
            p_label, p_acc, p_val = svm_predict(dummy_y, test_km, m, '-b 1 -q')
            pred_fp[:,i] = [p[0] for p in p_val]
        else:
            p_label, p_acc, p_val = svm_predict(dummy_y, test_km, m, '-b 0 -q')
            pred_fp[:,i] = p_label
    return pred_fp
    
