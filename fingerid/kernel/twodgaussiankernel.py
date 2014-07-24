"""
=============================================
Compute probability product kernel for MS/MS.
=============================================
"""

from mskernel import Kernel
import numpy

class TwoDGaussianKernel(Kernel):

    def __init__(self, sm, si):
        """ 
        One can use cross validation to find best sm and si. 

        """        
        Kernel.__init__(self)
        self.sm = sm # variance for mass
        self.si = si # variance for intensity

    def compute_train_kernel(self, spectra):
        """ 
          This function computes the probability product kernel 
          with peaks and mloss feature.
          
        """
        #print "computing training kernel now ..."
        #sm = self._get_mass_error(spectra)
        sm = self.sm
        si = self.si

        cleaned_spectra = spectra
        n = len(cleaned_spectra)
        
        km_peaks = numpy.zeros((n,n))
        km_mloss = numpy.zeros((n,n))
        for i in range(n):
            for j in range(i,n):
                spec_i = cleaned_spectra[i]
                spec_j = cleaned_spectra[j]

                i_peaks = self._peaks_to_matrix(spec_i)
                j_peaks = self._peaks_to_matrix(spec_j)
                km_value_peaks = self._gaussprodmixture(i_peaks,j_peaks,sm,si)
                
                km_peaks[i,j] = km_value_peaks
                km_peaks[j,i] = km_peaks[i,j]

                i_mloss = self._mloss_to_matrix(spec_i, spec_i.precursor)
                j_mloss = self._mloss_to_matrix(spec_j, spec_j.precursor)
                km_value_mloss = self._gaussprodmixture(i_mloss,j_mloss,sm,si)

                km_mloss[i,j] = km_value_mloss
                km_mloss[j,i] = km_mloss[i,j]
                #print "Computing kernel value for (%d, %d)" % (i,j)

        km_peaks = self._normalize_km(km_peaks)
        km_mloss = self._normalize_km(km_mloss)
        km = km_peaks + km_mloss

        km = self._normalize_km(km)
        return km

    def compute_test_kernel(self, test_spectra, train_spectra):
        """ 
          This function computes the probability product kernel 
          with peaks and mloss feature.
        """
        #print "computing test kernel now ..."

        sm = self.sm
        si = self.si


        cleaned_train_spectra = train_spectra
        cleaned_test_spectra = test_spectra

        n_train = len(cleaned_train_spectra)
        n_test = len(cleaned_test_spectra)
        km = numpy.zeros((n_test,n_train))

        for i in range(n_test):
            for j in range(n_train):
                #print "Computing kernel value for (%d, %d)" % (i,j)
                spec_i = cleaned_test_spectra[i]
                spec_j = cleaned_train_spectra[j]
                i_peaks = self._peaks_to_matrix(spec_i)
                j_peaks = self._peaks_to_matrix(spec_j)
                km_value_peaks = self._gaussprodmixture(i_peaks,j_peaks,sm,si)

                i_mloss = self._mloss_to_matrix(spec_i, spec_i.precursor)
                j_mloss = self._mloss_to_matrix(spec_j, spec_j.precursor)
                km_value_mloss = self._gaussprodmixture(i_mloss,j_mloss,sm,si)
                
                km_peaks_ii = self._gaussprodmixture(i_peaks,i_peaks,sm,si)
                km_mloss_ii = self._gaussprodmixture(i_mloss,i_mloss,sm,si)

                km_peaks_jj = self._gaussprodmixture(j_peaks,j_peaks,sm,si)
                km_mloss_jj = self._gaussprodmixture(j_mloss,j_mloss,sm,si)

                km_peaks_ij = km_value_peaks / numpy.sqrt(km_peaks_ii*km_peaks_jj)
                km_mloss_ij = km_value_mloss / numpy.sqrt(km_mloss_ii*km_mloss_jj)
                # normalize km
                km[i,j] =  (km_peaks_ij + km_mloss_ij) / 2
        return km


#    def _clean(self, spectra):
#        """  Normalize intensity and remove small peaks
#        Global normalization gives the same result as local normalization
#        which is also simpler, so use local normalization!
#        """
#        for spec in spectra:
#            sum_inten = 0
#            peaks_cleaned = []
#            for mass, inten in spec.peaks:
#                sum_inten = sum_inten + inten
#            for mass, inten in spec.peaks:
#                peaks_cleaned.append((mass,float(inten)/sum_inten))
#                spec.peaks = peaks_cleaned
#        return spectra
        
    def _normalize_km(self, km):
        n = len(km)
        for i in range(n):
            if km[i,i] == 0:
                km[i,i] = 1.0/100000
        return km / numpy.array(numpy.sqrt(numpy.mat(numpy.diag(km)).T * numpy.mat(numpy.diag(km))))  
        
#    def _get_mass_error(self, spectra): 
#        """ return the variance of the mass error """
#        error = []
#        for spec in spectra:
#            pre = spec.precursor
#            exact_mass = spec.mass
#            if abs(pre - exact_mass) > 2:  # do not use it in estimation
#                continue
#            if spec.mode == "POSITIVE":
#                error.append(pre - 1.007825 - exact_mass)
#            else:
#                error.append(pre + 1.007825 - exact_mass)
#        self.sm = abs(numpy.mean(numpy.square(error)))
#        return self.sm

    def _peaks_to_matrix(self, spec):
        """ convert list of peaks into ndarray form """
        diff =  spec.precursor - spec.mass # align the peaks. remove the addcut
        n = len(spec.peaks)
        matrix = numpy.zeros((n,2))
        count = 0
        for mass,inten in spec.peaks:
            matrix[count,0] = mass - diff # remove the addcut
            matrix[count,1] = inten
            count = count +1
        return matrix

    def _mloss_to_matrix(self, spec, pre): 
        n = len(spec.peaks)
        matrix = numpy.zeros((n,2))
        count = 0
        for mass,inten in spec.peaks:
            matrix[count,0] = abs(pre-mass)
            matrix[count,1] = inten
            count = count +1
        matrix[matrix == 0] = 0
        return matrix


    def _gaussprodmixture(self, i_peaks, j_peaks, sm, si):
        X1 = i_peaks; X2 = j_peaks;
        N1 = numpy.size(X1,0); N2 = numpy.size(X2,0)
        if N1 == 0 or N2 == 0:
            raise Exception("[ERROR]:No peaks when computing the kernel.(try not clean the peaks)")
        constant = 1.0/(N1*N2)*0.25/(numpy.pi*numpy.sqrt(sm*si))
        mass_term = 1.0/sm * numpy.power(numpy.kron(X1[:,0].flatten(),numpy.ones(N2)) - numpy.kron(numpy.ones(N1),X2[:,0].flatten()),2)
        inte_term = 1.0/si * numpy.power(numpy.kron(X1[:,1].flatten(),numpy.ones(N2)) - numpy.kron(numpy.ones(N1),X2[:,1].flatten()),2)
        #return constant*sum(numpy.exp(-0.5*(mass_term + inte_term)))
        return constant*sum(numpy.exp(-0.25*(mass_term)))

