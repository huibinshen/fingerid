"""
Search the molecule database for identification
"""
import numpy
import pickle
import os

def search(test_spectra, pred_fps, ppm, cv_accs_f, db, out_f, fpset='part'):
    """
    Search the molecular database using the predicted fingerprints

    Parameters:
    ----------
    test_spectra: the list of spectrum instance of unkown molecule
    pred_fp: numpy 2d array or the file contains it.
             The prediction of fingerprints
    ppm: search mass window
    cv_accs_f: string, file stores cross validation accuracy in training.
    db: string, which database to search, now only accept kegg
    out_f, string, the file contains the identification result
    """

    if db.lower() == "kegg":
        # get the dir in which the script being run
        script_dir = os.path.dirname(os.path.abspath(__file__))
        resource_dir = script_dir + "/../data"

        mass_db = resource_dir + "/../data/kegg_mass"
        if fpset == 'part':
            fp_db = resource_dir + "/../data/kegg_fp.dict"            
        elif fpset == 'all':
            fp_db = resource_dir + "/../data/kegg_all_fp.dict"
        else:
            raise Exception("[ERROR] Fingerprint set either 'part' or 'all'")   
    else:
        raise Exception("[ERROR] Only support searching Kegg now!")
    
    if type(pred_fps) == str:
        pred_fps = numpy.loadtxt(pred_fps)
    
    n = len(test_spectra)
    # read mass_db in
    ID_ARRAY, MASS_ARRAY = _get_id_mass(mass_db)
    # read fp_db in
    f = open(fp_db,"rb")
    FP_DB = pickle.load(f)
    f.close()
    # read the cross validation accuracies in
    cv_accs = numpy.fromfile(cv_accs_f,sep=" ")
    n_fp = len(cv_accs)

    # get the guess mass
    guessed_mass = _guess_mass(test_spectra)
#    real_mass = []
#    for spec in test_spectra:
#        real_mass.append(spec.mass)
#    real_mass = numpy.array(real_mass)
#    tmp = numpy.zeros((n,3))
#    tmp[:,0] = guessed_mass
#    tmp[:,1] = real_mass
#    tmp[:,2] = numpy.abs(guessed_mass - real_mass)
#    numpy.savetxt('mass_diff.txt',tmp,fmt="%.4f")

    all_res = []   
    for i in range(n): # begin search for each query mass spectrum # 
        # get candidates after mass filtering
        tol = float(ppm)/1000000*guessed_mass[i]
        candis = numpy.array(ID_ARRAY[abs(guessed_mass[i]-MASS_ARRAY)<=tol])
        pred_fp = pred_fps[i,:]; n_candi = len(candis)        
        one_res = _sort_candis(pred_fp,candis,n_fp,n_candi,FP_DB,cv_accs)
        all_res.append(one_res)

    w = open(out_f,"w")
    i = 0
    for one_res in all_res:
        one_res_str = []
        for iden, score in one_res:
           one_res_str.append(str(iden)+" "+str(score))
        w.write(test_spectra[i].f_name+" "+";".join(one_res_str)+"\n")
        i = i+1
    w.close()
    print "search done!"


def _sort_candis(pred_fp,candis,n_fp,n_candi,FP_DB,train_acc):
    """
    Sort the candidates based on fingerprint matching
    """
    if n_candi == 0:
        return []
    if n_candi == 1:
        return [(candis[0],1)]
    tmp_candis = [] # some candis are only in mass db but not in fp_db
    for i in range(n_candi):
        candi = candis[i]
        if candi not in FP_DB: # if not in FP_DB, continue  
            continue
        else:
            tmp_candis.append(candi)
    if len(tmp_candis) == 0:
        return []
    
    n_candi = len(tmp_candis)
    candis = numpy.array(tmp_candis)
    temp_matrix = numpy.zeros((n_candi,n_fp))
    scores = numpy.zeros(n_candi)
    res = [] # the result list of (ids, scores)  
    #construct temp matrix of fps for candis   
    for i in range(n_candi):
        candi = candis[i]
        db_fp = FP_DB[candi]
        # 528 bit string conver to numpy 1-d array
        db_fp = numpy.array(map(int,[db_fp[j] for j in range(len(db_fp))]))
        temp_matrix[i,:] = db_fp

    train_acc[train_acc==1] = 0.99999

    pred_fp_rep = numpy.tile(pred_fp,(n_candi,1))
    train_acc_rep = numpy.tile(train_acc,(n_candi,1))
    temp_matrix = temp_matrix * 2 - 1 # fp in db are 0/1 s

    scores = numpy.sum(numpy.log(train_acc_rep*(pred_fp_rep==temp_matrix)+
                           (1-train_acc_rep)*(pred_fp_rep!=temp_matrix)),1)
    sorted_candis = candis[numpy.argsort(scores)]
    sorted_scores = numpy.sort(scores)

    # map scores to [0,1]
    sorted_scores = sorted_scores + abs(sorted_scores[0]) + 1
    max_score = float(max(sorted_scores))
    sorted_scores = (sorted_scores+max_score/n_candi/n_candi) / (max_score + max_score/n_candi)
    
    for i in range(len(sorted_scores)-1,-1,-1):
        res.append((sorted_candis[i],sorted_scores[i]))
    return res


def _get_id_mass(mass_db):
    """
    Parse mass_db file (f_string) into two numpy arrays
    """
    ids = []; mass = []
    f = open(mass_db); count = 0
    while True:
        line = f.readline()
        if not line:
            break
        words = line.split()
        k = words[0]; m = float(words[1])
        ids.append(k); mass.append(m)
        count = count +1
    ids = numpy.array(ids).reshape((count,1))
    mass = numpy.array(mass).reshape((count,1))
    return ids,mass

def _guess_mass(spectra):
    """
    Guess the exact mass of every spectum in spectra
    Using the equation: exact_mass = precursor +/- 1.00785
    """
    guess_mass = []
    for spec in spectra:
        guess_mass.append(spec.mass)
#        continue
#        g_mass = 0
#        if spec.mode == "POSITIVE":
#            g_mass = (spec.precursor-1.00785)
#        elif spec.mode == "NEGATIVE":
#            g_mass = (spec.precursor+1.00785)
#        else:
#            raise Exception("ERROR: spectrum ion mode wrong!")
#        guess_mass.append(g_mass)
    return numpy.array(guess_mass)
