"""
===========================================
MS/MS parser for metlin data in .msx format
===========================================
"""

import glob
import pickle
import os
import numpy
import xml.etree.ElementTree as ET

from spectrum import Spectrum
from parser import Parser
from util import sortbyfilenames

class MetlinParser(Parser):
    """
    A metlin msx file parser.
    The spectrum instance for this parser has a extra field metlin_id
    The intensity will be normalized by 100 and 
    peaks less than intensity 1 will be removed.
    """

    def parse_file(self, f_path, kegg_mass, kegg_inchi):
        """ 
        parse file into Spectrum instance 
        """
        
        return self._parse_metlin_file(f_path, kegg_mass, kegg_inchi)

    def parse_dir(self, dir_path="NULL"):
        """ parse directory of MS/MS data into list of Spectrum instance """
        spec_list = []
        #dir_path = os.path.abspath(dir_path)

        # if user's path is not having a "/"
        if dir_path[-1] != "/": 
            dir_path = dir_path + "/"

        # invoke parse file for every file in the dir_path directory
        files = glob.glob(dir_path+'*.msx')
        
        # get the dir in which the script being run
        script_dir = os.path.dirname(os.path.abspath(__file__))
        resource_dir = script_dir + "/../data"

        f = open("%s/kegg_inchi.dict" % resource_dir,"rb")
        kegg_inchi = pickle.load(f)
        f.close()

        f = open("%s/kegg_mass" % resource_dir,"r") 
        lines = f.read().split("\n")
        f.close()
        kegg_mass = {}
        for line in lines:
            if len(line) <1:
                break
            k_id = line.split()[0]; k_mass = float(line.split()[1])
            kegg_mass[k_id] = k_mass

        for f in files:
            specs = self.parse_file(f, kegg_mass, kegg_inchi)
            spec_list = spec_list + specs

        return sortbyfilenames(spec_list, files)

    def _parse_metlin_file(self, f_path, kegg_mass, inchi):
        tree = ET.parse(f_path)
        root = tree.getroot()

        eles = root.findall("./ExperimentInformations/Comment")
        for ele in eles:
            if ele.get('Id')=='kegg':
                kegg_id = ele.get('Value')
            if ele.get('Id') == 'Metlin-ID':
                metlin_id = ele.get('Value')
            if ele.get("Id") == 'cas':
                cas = ele.get('Value')

        eles = root.findall("./ExperimentInformations")
        for ele in eles: # should have only one element
            mass_diff = float(ele.attrib['ModificationMass'])
            c_name = ele.attrib['CompoundName']
            c_formula = ele.attrib['MolecularFormula']

        spectra = []; spectra_ms1 = []
        eles = root.findall("./Spectra/Spectrum")
        for ele in eles:
            if ele.get("MSLevel") == "1":
                spectra_ms1.append(ele)
            if ele.get("MSLevel") == "2":
                spectra.append(ele)
        if kegg_id not in kegg_mass:
            print "Ignore %s:kegg_mass doesn't have %s" % (f_path,kegg_id)
            return []

        if kegg_id not in inchi:
            print "Ignore %s:kegg_inchi doesn't have %s" % (f_path,kegg_id)
            return []
        
        mass = kegg_mass[kegg_id]
        inchi = inchi[kegg_id]
        spectra_list = []
        for spec in spectra:
            ce = int(spec.attrib['CollisionEnergy'])
            spectrum = Spectrum()
            spectrum.f_name = f_path
            spectrum.mass = float(mass)
            spectrum.precursor = mass + mass_diff
            spectrum.mode = "POSITIVE"
            spectrum.inchi = inchi
            spectrum.cas = cas
            spectrum.pubchem_sid = "NULL"
            spectrum.pubchem_cid = "NULL"
            spectrum.kegg_id = kegg_id
            spectrum.metlin_id = metlin_id
            spectrum.ce = ce
            peaks = spec.findall("Peak")
            _peaks = []
            for peak in peaks:
                _mass = float(peak.get("Mass")); _inten = float(peak.get("Intensity"))
                if _inten > 1:
                    _peaks.append((_mass,_inten/100))
            spectrum.peaks = _peaks
            spectra_list.append(spectrum)
        return spectra_list



