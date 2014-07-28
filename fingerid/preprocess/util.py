
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

