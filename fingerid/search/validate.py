
import numpy

def validate(spectra, res_f, out_f):
    f = open(res_f)
    data = f.read()
    f.close()

    n = len(spectra)
    # columns is rank and the total candidate size
    # -1 means not in the list
    # -2 means no list produced
    # init res to -1
    res = numpy.zeros((2,n)) -1 
    lines = data.split('\n')
    count = 0
    for i in range(n):
        line = lines[i]
        line = line[line.find(' ')+1:]
        if len(line) == 0:
            res[0,count] = -2
            res[1,count] = -2
            count = count +1
            continue
        words = line.split(';')

        rank = 1
        for word in words:
            _id, score = word.split(' ')
            if _id == spectra[count].kegg_id:
                res[0,count] = rank
                res[1,count] = len(words)
              
            rank = rank +1
        count = count +1
    numpy.savetxt(out_f, res, fmt="%d")

