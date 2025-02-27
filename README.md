# OED
NOTE that this software has been adapted, modified and improved for different projects and aplications:
For more comprehensive READMEs go to the subfolders for each project

Compile models from a manually defined set of equations to AMICI, this applies to any ODE system
Perform identification analysis and optimal experimental design (i.e. obtain a series of time dependent inupts for you model which improve parameter estimations)

Beta version 1.0 (march 24 2022)

This is a beta version of the ODE compilation and experimental design program 
we will prepare a pip version in the near future (August 2022). 

More function will be added over time.

You can download the code folder to your desktop, this contains 2 subfolders

-Dependencies 

-Main

within the dependencies script you will find the write up of
	Modelobject/model builder and compilation with AMICI (you can make use of subfunctions in the 
	Measurementobject and experiment generation
	ModelSolver where model object with conditions is simulated
	ParameterOptimization (fitting data)
	ExperimentalOptization (optimizing control inputs)
	
In the demonstration folder are 2:

-A model folder, where we define a model or a set of models
-A main folder where we call the different optimization and analysis scripts

The user will need to install the following packages

-Tellurium  https://pypi.org/project/python-libsbml/](https://pypi.org/project/tellurium/

-Roadrunner https://pypi.org/project/libroadrunner/

-Libsbml    https://pypi.org/project/python-libsbml/    

-AMICI      https://amici.readthedocs.io/en/latest/about.html

-pyDOE      https://pythonhosted.org/pyDOE/ 

-sobolseq   https://pypi.org/project/sobol-seq/ 

-importlib  https://pypi.org/project/importlib/


The next version of this script will include a pip install for all dependencies excluding AMICI

for the AMICI installation we refer the user to the AMICI website https://amici.readthedocs.io/en/latest/about.html

Note that the windows install requires path changes, similarly amici requires specific compilers to be present.
Note that for the compilation of the model to AMICI we reload the script, when we create the SBML file from the model the user defines we use tellurium and libsbml, to our knowlegde these have conflicting packages with AMICI thus we reload the system (with packages ordered differently) to compile the model after the SBML file has been created.  


