""" 
======================================================
The class of internal representation of a tandem MS/MS
======================================================
"""

class Spectrum:

    """
    Attributes
    ----------
    
    f_name : string, default: "NULL"
        The name of the file containing the MS/MS.

    precursor : float, default: 0
        The precursor of the this tandam MS/MS.

    peaks: list of tuples [(mass,int), ...], default: []
        The data structure storing all the peaks

    Optional attributes
    -------------------
    mode: string, optional default: "POSITIVE"
        The ion mode of the mass spectrometer

    mass: float, optional, default: 0
        The exact mass the measured molecule

    inchi: string, optional, default = "NULL"
        In the traing, inchi is required for generate fingerprints.
        In teh testing, inchi can be empty.

    cas: string, optional, default: "NULL"
        The cas number of the measured molecule

    metlin_id: string, optional, default: "NULL"
        The metlin id of the measured molecule

    pubchem_sid: string, optional, default: "NULL"
        The pubchem substance id of the measured molecule

    pubchem_cid: string, optional, default: "NULL"
        The pubchem compound id of the measured molecule 

    kegg_id: string, optional, default: "NULL"
        The kegg compound id of the measured molecule

    ce: int, optional, default: 0
    
    """

    def __init__(self, f_name="NULL", precursor=0, peaks=[], mode="POSITIVE", mass=0, 
                 inchi="NULL", cas="NULL", metlin_id='NULL',pubchem_sid="NULL",
                 pubchem_cid="NULL", kegg_id="NULL", ce=0):
        
        self.f_name = f_name
        self.mass = mass
        self.precursor = precursor
        self.mode = mode
        self.peaks = peaks
        self.inchi = inchi
        self.cas = cas
        self.metlin_id = metlin_id
        self.pubchem_sid = pubchem_sid
        self.pubchem_cid = pubchem_cid
        self.kegg_id = kegg_id
        self.ce = ce

    def __str__(self):
        return "f_name:%s\nmass:%s\nprecursor:%s\nmode:%s\npeaks:%s\ninchi:%s\ncas:%s\npubchem_sid:%s\npubchem_cid:%s\nkegg_id:%s\nmetlin_id:%s\ncollision_energy:%s" % (self.f_name, self.mass, self.precursor, self.mode, self.peaks, self.inchi, self.cas, self.pubchem_sid, self.pubchem_cid, self.kegg_id, self.metlin_id,self.ce)

