"""
============================================================
This is the base class compute the msss spectrometry kernel.
============================================================
"""
import numpy

class Kernel:

    def __init__(self):
        # the maximum intensity of peaks observed in training data
        self._max_inten = 0 

    def compute_train_kernel(self, spectra):
        """ compute the kernel for spectra,
            to be implemented by inheritant class.
        
        Parameters
        ----------
        spectra: a list of spectrum instance

        Returns
        ------
        A numpy 2-d array with the dim as n*n, 
        (n is the number of spectra instance)

        """

        print "Should not be called with this class"

    def compute_test_kernel(self, test_spectra, train_spectra):
        """ compute the test kernel for spectra,
            to be implemented by inheritant class.
        
        Parameters
        ----------
        test_spectra: a list of spectrum instance need to be tested
        train_spectra: a list of training spectrum instance

        Returns
        -------
        A numpy 2-d array with the dim as m*n, 
        (n is the number of training spectra instance)
        (m is the number of testing spectra instance)
        """
        
        
    def read_kernel(self, kernel_path):
        """ read the kernel for spectra
        
        Parameters
        ----------
        kernel_path: kernel has been computed, read the kernel in.

        Returns
        -------
        A numpy 2-d array with the dim as n*n, 
        (n is the number of spectra instance)

        """

        return numpy.loadtxt(kernel_path)

    def write_kernel(self, kernel, write_path):
        """ write the kernel in a file
        
        Parameters
        ----------
        kernel: The numpy 2-d array with the shape n * n

        write_path: the file to write the kernel in

        """
        
        numpy.savetxt(write_path, kernel)
        
        
        
