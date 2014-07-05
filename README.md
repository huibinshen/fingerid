FingerID 1.4
============

This is FingerID 1.4 release. This version will only focus on the fingerprints
prediction, without fingerprints generation and compound database retrieval.

This package utilize fragmentation tree, another view of MS/MS spectra, to 
improve fingerprints prediction.

The previous versions are hosted on sourceforge: http://sourceforge.net/projects/fingerid/


Changes
=======

- Fingerprints generation and compound database retrieval are delegated to user.
- Cleaning of the MS/MS spectra is delegated to user.
- Merge of the MS/MS spectra is delegated to user.
- Fragmentation tree information is added in the input by kernel and multiple
kernel learning.

Dependencies
============

- Python >= 2.7
- Numpy >= 1.4.0
- [LibSVM](http://cvxopt.org/install/index.html) >= 3.17 python interface
- [cvxopt](http://cvxopt.org/install/index.html), optional, only needed if using 'ALIGNF' to combine the kernels

Install
=======
- If you have root permission:

  ```python
  python setup.py install
  ```

- or if you do not have root permission:

  ```python
  python setup.py install --user
  ```

- or in your python script (preferred):

  ```python
  import sys
  sys.path.append("path_to_this_foler")
  ```

Instructions
============

To use the package, three steps (parse, kernel, predict) are needed sequentially.
Two examples are also provided in shen_ISMB2014.py and train_test.py.

Parse
-----

Parse MS/MS spectra to the internal representation.

- For the MS/MS data in the format as example dataset provided in the package, one can use the following: 

  ```python
  from fingerid.preprocess.msparser import MSParser
  # ms_folder is the folder for all the spectra.
  ms_list = msparser.parse_dir(ms_folder) 
  ```

- For the MS/MS data downloaded from MassBank:

  ```python  
  from fingerid.preprocess.massbankparser import MassBankParser
  mbparser = MassBankParser()
  # ms_folder is the folder for all the spectra.
  ms_list = mbparser.parse_dir(ms_folder)
  ```

- For the MS/MS data downloaded from Metlin (.msx format):

  ```python
  from fingerid.preprocess.metlinparser import MetlinParser
  mlparser = MetlinParser()
  # ms_folder is the folder for all the spectra.
  ms_list = mlparser.parse_dir(ms_folder)
  ```

- For the fragmentation tree in .dot format (fgtree_folder is the folder name for fragmentation tree data):

  ```python
  from fingerid.preprocess.fgtreeparser import FragTreeParser
  fgtreeparser = FragTreeParser()
  trees = fgtreeparser.parse_dir(fgtree_folder)
  ```

Kernel
------

Two types of kernel functions are provided. For the MS/MS data, "PPK" kernel is used:
  
  ```python
  from fingerid.preprocess.msparser import MSParser
  from fingerid.kernel.twodgaussiankernel import TwoDGaussianKernel
  train_ms_list = msparser.parse_dir(train_ms_folder)
   
  # Compute the PPK kernel with m/z variance sm and intensity variance si.
  # In practice, tune the sm and si by cross validation is important.
  kernel = TwoDGaussianKernel(sm, si)
  train_km = kernel.compute_train_kernel(train_ms_list)

  # When have test data, to compute test kernel
  test_ms_list = msparser.parse_dir(test_ms_folder)
  test_km = kernel.compute_test_kernel(test_ms_list,train_ms_list)
  ```

For fragmentation tree:

- parse fragmentation tree

  ```python
  from fingerid.preprocess.fgtreeparser import FragTreeParser
  fgtreeparser = FragTreeParser()
  train_trees = fgtreeparser.parse_dir(train_fgtree_folder)
  ```

- Compute training kernel

  ```python
  kernel = FragTreeKernel()
  # Kernel can be "NB","NI","LB","LC","LI","RLB","RLI","CPC","CP2","CPK","CSC"
  train_tree_km = kernel.compute_kernel(train_trees, "NB")
  ```  

- When have test data for fragmentation trees

  ```python
  test_trees = fgtreeparser.parse_dir(test_fgtree_folder)  
  n_train = len(train_trees)
  n_test = len(test_trees)
  n = n_train + n_test
  trees = train_trees + test_trees
  kernel = FragTreeKernel()
  # Kernel can be "NB","NI","LB","LC","LI","RLB","RLI","CPC","CP2","CPK","CSC"
  tree_km = kernel.compute_kernel(trees, "NB")
  train_tree_km = tree_km[0:n_train, 0:n_train]
  test_tree_km = tree_km[n_train:n, 0:n_train]
  ```

To combine the kernel using MKL (UNIMKL, ALIGN, ALIGNF):

  ```python
  # km_list is a list of kernel matrices (numpy 2d array).
  # output is fingerprint matrix (numpy 2d array).
  # The MKL algorithms can be 'UNIMKL', 'ALIGN' and 'ALIGNF'.
  # ckm is combined kernel and kw is the weights for the kernels.
  # The weights can be used to combine the test kernel.
  from fingerid.kernel.mkl import mkl
  ckm, kw = mkl(km_list, output, 'ALIGN')
  ```

Predict
----------

To perform cross validation on training data:

  ```python
  from fingerid.model.internalCV_mp import internalCV
  # kernel is the kernel matrix (numpy 2d array)
  # labels is fingerprint matrix (numpy 2d array).
  # n_folds is the number of folds used in the cross validation
  # pred_f is the file to be written for the prediction
  # select_c is a boolean variable specify whether to do C selection in SVM.
  n_folds = 5
  pred_f = "prediction.txt"
  internalCV(kernel, labels, n_folds, pred_f, select_c=False)
  ```
 
To perform cross validation on training data with multiple processes. This is
useful when you have many fingerprints (output) to train:

  ```python
  from fingerid.model.internalCV_mp import internalCV_mp
  # n_p is the number of processes to be used
  internalCV_mp(kernel, labels, n_folds, pred_f, select_c=False, n_p=8)
  ```

To train the model on all the data instead of doing cross validation:

  ```python
  from fingerid.model.trainSVM import trainModels
  # model_f is the file to store the trained models
  models = trainModels(kernel, labels, select_c=False, n_p)
  ```

To predict on the test data using trained models:

  ```python
  from fingerid.model.trainSVM import trainModels
  from fingerid.model.predSVM import predModels
  model_dir = "MODELS" # model_dir is the folder to put the trained models
  trainModels(train_kernel, labels, model_dir, select_c=False, n_p)
  preds = predModels(test_kernel, n_fp, model_dir) # n_fp is the number of fingerprints
  ```

References
==========
Huibin Shen, Kai D\"ahrkop, Sebastian B\"ocker and Juho Rousu: Metabolite Identification through Multiple Kernel Learning on Fragmentation Trees. In the proceedings of ISMB 2014, Bioinformatics 30(12), i157-i164 (2014). 

Huibin Shen, Niocola Zamboni, Markus Heinonen, Juho Rousu: Metabolite identification through machine learning -- tackling casmi challenge using fingerid. Metabolites 3(2), 484--505 (2013).

Markus Heinonen, Huibin Shen, Niocola Zamboni, Juho Rousu: Metabolite identification and molecular fingerprint prediction through machine learning. In the proceedings of MLSB 2012, Bioinformatics 28(18), 2333--2341 (2012).
