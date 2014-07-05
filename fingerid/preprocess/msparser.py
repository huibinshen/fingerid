"""
================================================
Parser of MS/MS in the format as example dataset
================================================
"""

import re
import commands
import os
import numpy

from spectrum import Spectrum
from parser import Parser

class MSParser(Parser):
    """
    The parser for ms format used in the example dataset.

    """

    def parse_file(self, f_path="NULL"):
        """ parse file into Spectrum instance """
        if f_path == "NULL":
            raise Exception("ERROR: please specify tandam MS/MS file path")
        return self._parse_ms_file(f_path)

    def parse_dir(self, dir_path="NULL"):
        """ parse directory of MS/MS data into list of Spectrum instance 
        """
        
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

    def _parse_ms_file(self, f_path):
#        print "Parse file:",f_path
        # read ms/ms file in                                                   
        f = open(f_path)
        data = f.read()
        f.close()

        # create Spectrum instance                                             
        spectrum = Spectrum(f_path)

        # set f_name                                                           
        spectrum.f_name = f_path
        # set metlin id
        spectrum.metlin_id = f_path[f_path.find("pos")+3:f_path.find(".")]
        # set precursor
        _precursor = re.findall("parentmass[: ]+([0-9\.]+)",data)
        if len(_precursor) > 0:
            precursor = float(_precursor[0])
        else:
            raise Exception("ERROR: precursor not set for %s!" % f_path)
        spectrum.precursor = precursor
        spectrum.mass = precursor - 1.00794
        # set peaks and intensity
        _peaks = []
        seg = False
        for line in data.split('\n'):
            if line.find("collision") != -1:
                seg = True
                continue
            if not line:
                seg = False
                continue
            if seg:
                words = line.split()
                mass = float(words[0])
                inten = float(words[1])
                _peaks.append((mass,inten))
        spectrum.peaks = _peaks

        return spectrum
