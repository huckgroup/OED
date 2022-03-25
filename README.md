# OED
compile models and perform identification analysis and optimal experimental design


Beta version 1.0 (march 24 2022)

This is a beta version of the ODE compilation and experimental design program 
we will prepare a pip version in the near future. 

We are building this pipeline to include ever more functions and ease of use

You can download the code folder to your desktop, this contains 2 subfolders

-Dependencies 
-Demonstration

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
-Tellurium
-Libsbml
-pyDOE
-sobolseq
-AMICI

The next version of this script will include a pip install for all dependencies excluding amici

for the AMICI installation we refer the user to the AMICI website https://amici.readthedocs.io/en/latest/about.html

Note that the windows install requires path changes, similarly amici requires specific compilers to be present.


