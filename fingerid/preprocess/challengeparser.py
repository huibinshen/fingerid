""" 
===================================
The parser for CASMI challenge 2012
===================================

"""
import re
import commands
import os 
import numpy
 
from spectrum import Spectrum
from parser import Parser

class ChallengeParser(Parser):

    """
    A massbank text file parser.
    """
    def parse_file(self, f_path="NULL"):
        """ parse file into Spectrum instance """
        if f_path == "NULL":
            raise Exception("ERROR: please specify tandam MS/MS file path")
        return self._parse_massbank_file(f_path)

    def parse_dir(self, dir_path="NULL"):
        """ parse directory of MS/MS data into list of Spectrum instance """
        spec_list = []
        dir_path = os.path.abspath(dir_path)

        # if user's path is not having a "/"
        if dir_path[-1] != "/": 
            dir_path = dir_path + "/"

        # invoke parse file for every file in the dir_path directory
        files = commands.getoutput("ls %s" % dir_path).split()
        for f in files:
            spec = self.parse_file(dir_path + f)
            spec_list.append(spec)

        return spec_list

    def _parse_massbank_file(self, f_path):
        print "Parse file:",f_path
        # read ms/ms file in
        f = open(f_path)
        data = f.read()
        f.close()
        
        # create Spectrum instance 
        spectrum = Spectrum(f_path)

        # set f_name
        spectrum.f_name = f_path

        # set precursor
        _precursor = re.findall("MS\$FOCUSED_ION: PRECURSOR_M/Z[: ]+([0-9\.]+)",data) 
        if len(_precursor) > 0:
            precursor = float(_precursor[0])
        else:
            _basepeak = re.findall("MS\$FOCUSED_ION: BASE_PEAK[: ]+([0-9\.]+)",data)
            if len(_basepeak)>0 :
                print ("WARNING: using base peak as precursor for %s!" % f_path)
                precursor = float(_basepeak[0])
            else:
                raise Exception("ERROR: precursor not set for %s!" % f_path)
        spectrum.precursor = precursor

        # set ion mode
        _mode = re.findall("ION_MODE ([A-Z]+)", data)
        if len(_mode) > 0:
            mode = _mode[0]
        else:
            _mode = re.findall("MODE ([A-Z]+)", data)
            if len(_mode)>0:
                print ("WARNING: ion mode is set by MODE for %s!" % f_path)
                mode = _mode[0]
            else:
                raise Exception("ERROR: mode not set for %s!" % f_path)
        spectrum.mode = mode

        if spectrum.mode == 'POSITIVE':
            spectrum.mass = spectrum.precursor - 1.00794
        else:
            spectrum.mass = spectrum.precursor + 1.00794

        _ppm = re.findall("SE\$SEARCH_PPM[: ]+([0-9\.]+)",data) 
        if len(_ppm) > 0:
            ppm = int(_ppm[0])
        else:
            raise Exception("ERROR: PPM not set for %s!" % f_path)            
        spectrum.ppm = ppm

        # set peaks
        _peaks = []
        lines = data.split("\n"); ready = False
        for line in lines:
            if len(line) == 0:
                continue
            if line.find("PK$PEAK") != -1:
                ready = True
                continue
            if ready:
                if line.find("N/A") != -1:
                    raise Exception("ERROR: no peaks in %s" % f_path)
                words = line.split()
                mass = float(words[0])
                inten = float(words[1])
                #mass = mass+numpy.random.normal(0,1e-8,1) # add noise
                #mass = float("%.3f" % mass)
                _peaks.append((mass,inten))
        spectrum.peaks = _peaks


        _ce = re.findall("COLLISION_ENERGY (\w+)",data)
        if len(_ce) > 0:
            ce = _ce[0]
            ce = ce.replace("eV","")
            if ce.isdigit():
                spectrum.ce = int(ce)
        return spectrum


    # http://stackoverflow.com/q/3844948/, check whether all the entry in lst are same
    #def _checkEqualIvo(self, lst):
    #    return not lst or lst.count(lst[0]) == len(lst)

