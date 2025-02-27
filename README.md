# OED
Optimal experimental design software for biochemical systems:
NOTE that this software has been adapted, modified and improved for different projects and aplications:
FOR MORE COMPREHENSIVE READMES GO TO THE SPECIFIC SUBFOLDERS

Compile models from a manually defined set of equations to AMICI, this applies to any ODE system
Perform identification analysis and optimal experimental design (i.e. obtain a series of time dependent inupts for you model which improve parameter estimations).

The user will need to install the following packages

-Tellurium  https://pypi.org/project/python-libsbml/](https://pypi.org/project/tellurium/

-Roadrunner https://pypi.org/project/libroadrunner/

-Libsbml    https://pypi.org/project/python-libsbml/    

-AMICI      https://amici.readthedocs.io/en/latest/about.html

-pyDOE      https://pythonhosted.org/pyDOE/ 

-sobolseq   https://pypi.org/project/sobol-seq/ 

-importlib  https://pypi.org/project/importlib/

for the AMICI installation we refer the user to the AMICI website https://amici.readthedocs.io/en/latest/about.html

Note that the windows install requires path changes, similarly amici requires specific compilers to be present.
Note that for the compilation of the model to AMICI we reload the script, when we create the SBML file from the model the user defines we use tellurium and libsbml, to our knowlegde these have conflicting packages with AMICI thus we reload the system (with packages ordered differently) to compile the model after the SBML file has been created.  


