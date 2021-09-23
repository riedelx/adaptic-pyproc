# ADAPTIC PyProc

*ADAPTIC PyProc* is a repository that contains useful modules written in Python for pre-processing and post-processing data for finite element analysis software [ADAPTIC](https://www.imperial.ac.uk/media/imperial-college/research-centres-and-groups/computational-structural-mechanics/ADAPTIC_Manual.pdf). 

## Introduction

[ADAPTIC](https://www.imperial.ac.uk/media/imperial-college/research-centres-and-groups/computational-structural-mechanics/ADAPTIC_Manual.pdf) is a finite element analysis (FEA) software developed by Prof. Izzuddin at Imperial College London. It is an adaptive static and dynamic analysis program, which accounts for material and geometrical nonlinearities. It has been successfully used in the analyses of extreme loadings such as earthquakes, blast explosions, fire and sudden column removals. Several special structural forms have been developed including RC and composite slabs, cables, membranes and curved shells. Some of its novel features and findings using the software have been published in the leading international scientific journals (see [publications](http://imperial.ac.uk/people/b.izzuddin/publications.html)).

There is no modelling GUI available and it requires specially prepared input files. Once the analysis is completed, output files are produced and two graphics post-processing applications are available:

* _ADAPTIC_graphs_ for plotting X-Y graphs.
* _ADAPTIC_shapes_ for plotting deflected shapes.

The data syntax might be daunting for the new users, whilst for the experienced ones it could be time-consuming. *ADAPTIC PyProc* or aims to make the modelling and analysis easier. Processing of 2D and 3D models is available.


## Navigation and modules arrangement

The repository is divided into four sub-folders:
1. assets
2. libraries
3. documentation
4. models

### 1. Assets

This folder contains images and miscellaneous data used to obtain certain properties of the model.

### 2. Libraries 

These are all the python libraries _.py_ that provide the necessary classes and functions to run the analyses. These were designed for pre-processing and post-processing data in _ADAPTIC_. 

1. _utils.py_ - contains miscellaneous functions and classes used through the entire project.
2. _materials.py_ - material classes. Currently the following material models are available: stl1, con1, con2, stmdl2 and bond.
3. _sections.py_ - section classes. Currently the following section models are available: rss, rccs, rcts and isec.
6. _adaptic.py_ - this defines the *adaptic* class for post-processing.

### 3. Documentation

The documentation of the libraries is covered in the corresponding _.ipynb_ notebooks. 

1. _materials.ipynb_ - covers the definition and application of different material classes available in _materials.py_.
2. _sections.ipynb_ - describes different section classes available in _sections.py_. 

### 4. Models

Sample models were created to show how *ADAPTIC PyProc* could be incorporated into a project. The definition of a model is done in _preprocessor.ipynb_. When running this notebook a _.dat_ file is created and used as an input for _ADAPTIC_ (see [software documentation](https://spiral.imperial.ac.uk/handle/10044/1/4228)). The output files have to be then placed in the model folder to run postprocessor.ipynb_ notebook for post-processing. Each model contains the following folders and files:

2. _preprocessor.ipynb_ - this notebook generates the _.dat_ and _pkl_ files. It reads the libraries defined in _modules_ folder and modifies certain properties if necessary. It does not provide the overview of the model, which is covered in _analysis.ipynb_ instead.
2. postprocessor.ipynb_ - this notebook describes the model, unpickles the objects from pre-processing stage and reads the *ADAPTIC* output from the _num_ file. Next, post-processing of the data is carried out with an interpretation of the results.
5. _.dat_ file - input file for *ADAPTIC* analysis.
6. _.num_ file - output file from *ADAPTIC* that contains the numerical results at all requested load/time steps.
7. _.log_ file - output file from _ADAPTIC_ that contains solution progress log (not required by *ADAPTIC PyProc*).
8. _.plt_ file - *ADAPTIC* output plot file used by the native ADAPTIC post-processing programs (not required by *ADAPTIC PyProc*).

## Units

The units required in program are shown in the table below, unless stated otherwise.

| Unit  |
| :---: |
|  mm   |
|   N   |
|  MPa  |
| N/mm  |
| N/mm³ |

