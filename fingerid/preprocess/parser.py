"""
=======================================================================
The base class for parser.
=======================================================================

All the parsers which parse the tandam MS/MS spectra should implemented 
the parse_file and parse_dir function. 

"""

class Parser:

    def parse_file(self, f_path="NULL"):
        """ Parse the file in the parameter into the Spectrum instance 

        Parameters
        ----------
        f_path: string, the file location

        Return
        ------
        A instance of Spectrum
        """

        print "Should not be called within this class"

    def parse_dir(self, dir_path="NULL"):
        """ Parse the dir in the parameter into the list of Spectrum instance

        Parameters
        ----------        
        dir_path: string, the directory location

        Returns
        -------
        A list of instance of Spectrum
        """
        print "Should not be called within this class"

#    def _merge_spectra(self, spec_list):
#        """ 
#        Take all the peaks with same cas number (same molecule) in one peak list.
#        If the spectra does not have cas number, merge with kegg_id.
#        If the spectra does not have anything metioned above, just put it into the final list
        
#        Parameters
#        ----------
#        spec_list:, a list of spectrum instance
#        
#        Return
#        ------
#        A list of merged spectrums
#        """

#        print "STATUS: Merging spectra begin ..."
#        print "ATTENTIION: Spectra order will be different with file system after merging!"
#        cas_spectra = {}
#        merged_list = []
#        waiting_spectra = []
#        # for the same molecule, put all spectra (Spectrum instance) in the same list
#        merge_n = []
#        for spec in spec_list:
#            if spec.cas == "NULL":
#                waiting_spectra.append(spec)
#                continue
#            if spec.cas not in cas_spectra:
#                cas_spectra[spec.cas] = []
#            cas_spectra[spec.cas].append(spec)

#        for cas, cas_specs in cas_spectra.items():
#            # new spectrum has the same info except peaks and precursor
#            new_spec = cas_specs[0] 
#            peaks = []
#            prec = 0
            # take precursor closest to the molecule mass as precursor of new merged spectrum
#            prec_mass_diff = 1500 # for metabolites, 1500 is enough
#            for spec in cas_specs:
#                peaks = peaks + spec.peaks
#                diff = abs(spec.precursor - spec.mass)
#                #print spec.kegg_id, 'mass:', spec.mass, 'pre:',spec.precursor,spec.f_name
#                if diff < prec_mass_diff:
#                    prec_mass_diff = diff
#                    prec = spec.precursor
#                #print 'pre:',prec
#                #print
#            merge_n.append(len(cas_specs))
#            new_spec.peaks = peaks
#            new_spec.precursor = prec
#            merged_list.append(new_spec)

#        if len(waiting_spectra) < 1: # all the spectra have been merged
#            print "STATUS: Merging spectra done! (%d molecules, merge rate %f)" % (len(merged_list), float(sum(merge_n))/len(merge_n))
#            return merged_list
            
#        # for those spectra without cas try merging with kegg_id 
#        kegg_spectra = {}
#        new_waiting_list = []
#        for spec in waiting_spectra:
#            if spec.kegg_id == "NULL":
#                new_waiting_list.append(spec)
#                continue
#            if spec.kegg_id not in kegg_spectra:
#                kegg_spectra[spec.kegg_id] = []
#            kegg_spectra[spec.kegg_id].append(spec)
#        for kegg_id, kegg_specs in kegg_spectra.items():
#            # modes of spectra for the same molecule aren't the same, stop merging
#            #mode_list = [spec.mode for spec in cas_specs]
#            #if not self._checkEqualIvo(mode_list):
#            #    print("WARNING: mode different, don't merge and proceed")
#            #    return spec_list
#            
#            # new spectrum has the same info except peaks and precursor
#            new_spec = kegg_specs[0] 
#            peaks = []
#            prec = 0
#            # take precursor closest to the molecule mass as precursor of new spectrum
#            prec_mass_diff = 1000 # for metabolites, 1000 is enough
#            for spec in kegg_specs:
#                peaks = peaks + spec.peaks
#                diff = abs(spec.precursor - spec.mass)
#                if diff < prec_mass_diff:
#                    prec_mass_diff = diff
#                    prec = spec.precursor
#            new_spec.peaks = peaks
#            new_spec.precursor = prec
#            merged_list.append(new_spec)
#            merge_n.append(len(kegg_specs))
#        if len(new_waiting_list) < 1:
#            print "STATUS: Merging spectra done! (%d molecules, merge rate %f)" % (len(merged_list), float(sum(merge_n))/len(merge_n))
#            return merged_list        

        #for spec in new_waiting_list:
        #    print "%s has no cas and kegg_id, remove it!" % spec.f_name
        #    #merged_list.append(spec)
#        print "%s spectra have no cas and kegg_id" % len(new_waiting_list)

#        print "STATUS: Merging spectra done! (%d molecules, merge rate %f)" % (len(merged_list), float(sum(merge_n))/len(merge_n))
#        
#        return merged_list

