import numpy 

def writeIDs(fname, objects):
    """
    Write objects id into the file named by fnam

    Parameters:
    -----------
    fname: string, name of the file to be written 
    objects: list, list of objects (spectrum instances or fgtree instances)

    """
    files = []
    for obj in objects:
        files.append(obj.f_name)
    w = open(fname,"w")
    w.write("\n".join(files))
    w.close()

def sortbyfilenames(spectra, fnames):
    """ sort the spectrum in spectra by ids """
    return [spectra[i] for i in sorted(range(len(fnames)), key=lambda k: fnames[k])]

def centerTestKernel(km):
    """ centering test kernel when test kernel is not symmetric """
    nx, ny = km.shape
    ckm = numpy.zeros((nx,ny))
    for i in range(nx):
        for j in range(ny):
            ckm[i,j] = km[i,j] - numpy.mean(km[:,j]) - numpy.mean(km[i,:]) + numpy.mean(km)
    return ckm

