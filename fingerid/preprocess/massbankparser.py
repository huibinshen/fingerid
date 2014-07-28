""" 
==============================================
The parser parse MassBank format MS/MS spectra
==============================================

"""
import re
import commands
import os 
import numpy
 
from spectrum import Spectrum
from parser import Parser
from util import sortbyfilenames

class MassBankParser(Parser):
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
        fnames = []
        for f in files:
            spec = self.parse_file(dir_path + f)
            fnames.append(spec.fname)
            spec_list.append(spec)
        return sortbyfilenames(spec_list,fnames)

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

        # set mass
        _mass = re.findall("CH\$EXACT_MASS[: ]+([0-9\.]+)", data)
        if len(_mass) > 0:
            mass = float(_mass[0])
        else:
            raise Exception("ERROR: mass filed error in file %s " % f_path)
        spectrum.mass = mass

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

        # set inchi
        _inchi = re.findall("IUPAC: (.+)",data)
        if len(_inchi) > 0:
            if _inchi[0].find('unknown') != -1:
                print f_path, 'has no inchi!'
                inchi = _inchi[0]
                #raise Exception("Error: no inchi for %s!" % f_path)
            else:
                inchi = _inchi[0]
        else:
            raise Exception("Error: no inchi for %s!" % f_path)
        if "InChI=" not in inchi: # some inchi may not contains the head 
            inchi = "InChI=" + inchi
        spectrum.inchi = inchi

        # below are optional field for Spectrum
        _cas = re.findall("CH\$LINK: CAS[: ]+([0-9\-]+)", data)
        if len(_cas) > 0:
            cas = _cas[0]
            spectrum.cas = cas

        _metlin = re.findall("CH\$LINK: METLIN[: ]+([0-9]+)", data)
        if len(_metlin) > 0:
            metlin = _metlin[0]
            spectrum.metlin_id = metlin
        else:
            spectrum.metlin_id = 'NULL'

        _sid = re.findall("PUBCHEM SID[: ]+(\w+)", data)
        if len(_sid) > 0:
            sid = _sid[0]
            spectrum.pubchem_sid = sid
        else:
            _sid = re.findall("PUBCHEM[: ]+([0-9]+)", data)
            if len(_sid) > 0:
                sid = _sid[0]
                spectrum.pubchem_sid = sid

        _cid = re.findall("PUBCHEM CID[: ]+(\w+)", data)
        if len(_cid) > 0:
            cid = _cid[0]
            spectrum.pubchem_cid = cid

        _kegg_id = re.findall("LINK: KEGG (\w+)", data)
        if len(_kegg_id) > 0:
            kegg_id = _kegg_id[0]
            spectrum.kegg_id = kegg_id

        _ce = re.findall("COLLISION_ENERGY (\w+)",data)
        if len(_ce) > 0:
            ce = _ce[0]
            ce = ce.replace("eV","")
            if ce.isdigit():
                spectrum.ce = int(ce)
        return spectrum



