"""
===========================================================
Find the combined kernel matrix using UNIMKL, ALIGN, ALIGNF
===========================================================
"""
import numpy        

def mkl(km_list, labels, c_type):
    """Multiple kernel learning with UNIMKL, ALIGN, ALIGNF.

    The learned weights for kernels and the combined kernel are written out.
    
    Parameters:
    -----------
    km_list: list of numpy 2d array, a list of kernel matrix

    labels: numpy 2d array, the output

    c_type: str, the type of MKL, must be one of'UNIMKL', 'ALIGN', 'ALIGNF'

    Returns:
    --------

    train_km, numpy 2d array, combined kernel

    """
    print "Computing combined kernel for", c_type
    n_km = len(km_list)

    if c_type == 'UNIMKL':
        train_km = numpy.zeros(km_list[0].shape)
        for km in km_list:
            train_km = train_km + km
        train_km = normalize_km(train_km)
        w = np.ones(n_km)/n_km
        #numpy.savetxt(out_f, train_km)

    elif c_type == 'ALIGN':
        """Two stage model, weight as idenpedently centered alignment score"""
        ky = numpy.dot(labels,labels.T)
        ky_c = center(ky)
        ky_c = normalize_km(ky_c)
        ky_n = f_norm(ky_c)
        w = []
        train_km = numpy.zeros(km_list[0].shape)
        for km in km_list:
            km_c = center(km)
            km_c = normalize_km(km_c)
            s = f_dot(km_c, ky_c) / f_norm(km) / ky_n
            w.append(s)
            train_km = train_km + s*km_c
        train_km = normalize_km(train_km)
        w = numpy.array(w)
        #numpy.savetxt('kernel_weights_%s.txt' % c_type, numpy.array(w))
        #numpy.savetxt(out_f, train_km)

    elif c_type == 'ALIGNF':
        ky = numpy.dot(labels,labels.T)
        ky_c = center(ky)
        ky_c = normalize_km(ky_c)
        a = []
        kmc_list = []
        for km in km_list:
            km_c = center(km)
            km_c = normalize_km(km_c)
            a.append(f_dot(km_c, ky_c))
            kmc_list.append(km_c)
        a = numpy.array(a,dtype='d')
        n_km = len(kmc_list)
        M = numpy.zeros((n_km, n_km),dtype='d')
        for i in range(n_km):
            for j in range(i,n_km):
                M[i,j] = f_dot(kmc_list[i], kmc_list[j])
                M[j,i] = M[i,j]
        try:
            from cvxopt import matrix
            from cvxopt import solvers
        except ImportError:
            raise ImportError("Module <cvxopt> needed for ALIGNF.")

        # define variable for cvxopt.qp
        P = matrix(2*M)
        q = matrix(-2*a)
        G = matrix(numpy.diag([-1.0]*n_km))
        h = matrix(numpy.zeros(n_km,dtype='d'))
        sol = solvers.qp(P,q,G,h)
        w = sol['x']
        w = w/numpy.linalg.norm(w)
        # compute train_km with best w
        train_km = numpy.zeros(km_list[0].shape)
        for i in range(n_km):
            train_km = train_km + w[i]*kmc_list[i]
        train_km = normalize_km(train_km)
        #numpy.savetxt('kernel_weights_%s.txt' % c_type, numpy.array(w))
        #numpy.savetxt(out_f, train_km)
    else:
        raise Exception("MKL only supports UNIMKL, ALIGN, ALIGNF")
    return train_km, w

def center(km):
    """ centering km """
    m = len(km)
    I = numpy.eye(m)
    one = numpy.ones((m,1))
    t = I - numpy.dot(one,one.T)/m
    return numpy.dot(numpy.dot(t,km),t)

def f_norm(k):
    """ Compute frobenius norm of the matrix k """
    return numpy.sqrt(numpy.trace(numpy.dot(k.T,k)))

def f_dot(k1,k2):
    """ Compute frobenius product of matrix k1, k2 """
    return numpy.trace(numpy.dot(k1.T,k2))

def normalize_km(km):
    n = len(km)
    for i in range(n):
        if km[i,i] == 0:
            km[i,i] = 1.0/100000
    return km / numpy.array(numpy.sqrt(numpy.mat(numpy.diag(km)).T * numpy.mat(numpy.diag(km))))

