"""
================================================
Generate fingperprints for inchi using OpenBabel
================================================



"""
import commands
import numpy
import re
import os
import sys

class FPGenerator:
    """
    Attibutes
    ---------
    fp_type : string, default: "OpenBabel"
    fpset: string, either "part" or "all". "part" will generate FP3, FP4, MACCS,
    all together 528 bits while "all" will generate FP2, FP3, FP4, MACCS all 
    together 1856 bits.

    Functions:
    ----------
    generateFP(spectra: spectrum instance)

    """

    def __init__(self, fp_type="OpenBabel",fpset="part"):
        """

        Parameters:
        -----------
        fp_type, string, optional, current only support "OpenBabel".
        fpset: string, optional, generate 528 bits or full 1856 bits.

        """
        self.fp_type = fp_type
        self.fp_set = fpset
        if fp_type == "OpenBabel":
            if not os.environ.has_key('BABELPATH'):
                raise Exception("ERROR: BABELPATH are needed!") 
            else:
                self.babel_path =os.environ['BABELPATH'][1:]
        else:
            raise Exception("ERROR: Only allow OpenBabel fingerprints")

    def generateFP(self, spectra):
        """
        Parameter
        --------
        spectra: a list of spectra needs fingerprints
        fpset: 'part' (FP3,FP4,MACCS 528 bits) or 'all' (include FP2, 1856 bits)

        Returns
        -------
        fp: numpy.2d array [n_spectra,n_fp]

        """
        print "STATUS: Generating fingerprints now ..."
        if self.fp_type == "OpenBabel":
            return self._openbabel_fp_generator(self.babel_path, 
                                                spectra, self.fp_set)
        else:
            raise Exception("ERROR: Only allow OpenBabel fingerprints")
            
    def read_fp(self, file_path):
        """ read the fingerprint matrix in """
        return numpy.loadtxt(file_path)

    def write_fp(self, fp, write_path):
        """ write the fingerprint matrix in the file"""
        numpy.savetxt(write_path, fp, fmt="%d")
        
    def _openbabel_fp_generator(self, babel_path, spectra, fpset):
        n = len(spectra)
        if fpset == 'part':
            output = numpy.zeros((n,528))
            for i in range(n):
                s = spectra[i]
                inchi = s.inchi
                binstr = self._gen_528_babel_fp(babel_path, inchi)
                output[i,] = map(int, " ".join(binstr).replace("0","-1").split(" "))
        elif fpset == 'all':
            output = numpy.zeros((n,1856))
            for i in range(n):
                s = spectra[i]
                inchi = s.inchi
                binstr = self._gen_all_babel_fp(babel_path, inchi)
                output[i,] = map(int, " ".join(binstr).replace("0","-1").split(" ")) 
        else:
            raise('[ERROR]: can only generate openbale fingerprints!')
        return output

    def _gen_all_babel_fp(self, babel_path, inchi):
        mid = "temp"
        infn = mid+".inchi"; ofile = mid+".fpt"
        f = open(infn,"w")
        f.write(inchi)
        f.close()

        res = commands.getoutput("%s %s -ofpt -xfFP2 -h > %s" % (
                           babel_path+"/babel", infn, ofile))      
        commands.getoutput("%s %s -ofpt -xfFP3 -h >> %s" % (
                           babel_path+"/babel", infn, ofile))      
        commands.getoutput("%s %s -ofpt -xfFP4 -h >> %s" % (
                           babel_path+"/babel", infn, ofile))
        commands.getoutput("%s %s -ofpt -xfMACCS -h >> %s" % (
                           babel_path+"/babel", infn, ofile))
        #fp3 fp4 maccs has 55, 307, 166 functional groups        
        #featcounts = [1020, 55, 307, 166] 
        f = open(ofile)
        lines = f.readlines()
        lines = map(str.strip, lines)
        f.close()
        if len(lines) < 1:
            return "NULL"
        # the file contains three blocks on fingerprints:       
        # FP3, FP4 and MACCS       
        checksums = map(int, re.findall("([0-9]+) bits", " ".join(lines)))
        data = filter(lambda x: not x.startswith(">"), lines)
        fp2 = data[0:6]
        fp3 = data[6]
        fp4 = data[7:10]
        maccs = data[10:12]
        binstr = ["","","",""]
        # fp2
        hexstr = " ".join(fp2)
        for block in hexstr.split():
            binstr[0] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[0] = binstr[0][-1024:]
        # fp3
        hexstr = fp3
        for block in hexstr.split():
            binstr[0] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[1] = binstr[1][-55:]

        hexstr = " ".join(fp4)
        for block in hexstr.split():
            binstr[1] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[2] = binstr[2][-307:]

        hexstr = " ".join(maccs)
        for block in hexstr.split():
            binstr[2] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[3] = binstr[3][-166:]
        binstr = binstr[3]+ binstr[2] + binstr[1] + binstr[0]

        binstr = binstr[::-1]
        assert sum(checksums) == binstr.count("1")
        commands.getoutput("rm -f %s %s" % (infn, ofile))
        return binstr

    def _gen_528_babel_fp(self, babel_path, inchi):
        mid = "temp"
        infn = mid+".inchi"; ofile = mid+".fpt"
        f = open(infn,"w")
        f.write(inchi)
        f.close()
        commands.getoutput("%s %s -ofpt -xfFP3 -h > %s" % (
                           babel_path+"/babel", infn, ofile))      
        commands.getoutput("%s %s -ofpt -xfFP4 -h >> %s" % (
                           babel_path+"/babel", infn, ofile))
        commands.getoutput("%s %s -ofpt -xfMACCS -h >> %s" % (
                           babel_path+"/babel", infn, ofile))
        #fp3 fp4 maccs has 55, 307, 166 functional groups        

        featcounts = [55, 307, 166] 
        f = open(ofile)
        lines = f.readlines()
        lines = map(str.strip, lines)
        f.close()
        # the file contains three blocks on fingerprints:       
        # FP3, FP4 and MACCS       
        checksums = map(int, re.findall("([0-9]+) bits", " ".join(lines)))
        data = filter(lambda x: not x.startswith(">"), lines)

        fp3 = data[0]
        fp4 = data[1:4]
        maccs = data[4:6]
        binstr = ["","",""]

        # fp3
        hexstr = fp3
        for block in hexstr.split():
            binstr[0] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[0] = binstr[0][-55:]
        

        hexstr = " ".join(fp4)
        for block in hexstr.split():
            binstr[1] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[1] = binstr[1][-307:]

        hexstr = " ".join(maccs)
        for block in hexstr.split():
            binstr[2] += bin(int(block, 16))[2:].rjust(32, "0")
        binstr[2] = binstr[2][-166:]

        binstr = binstr[0] + binstr[1] + binstr[2]
#        binstr = binstr[2] + binstr[1] + binstr[0]
#        binstr = binstr[::-1]
        assert sum(checksums) == binstr.count("1")
        commands.getoutput("rm -f %s %s" % (infn, ofile))
        return binstr
