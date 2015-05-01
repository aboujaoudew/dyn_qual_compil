# dyn_qual_compil

Author: Wassim Abou-Jaoudé
Date:30 April 2015
------------------------------------------------------------------------------------------

Matlab program for the automatic generation of dynamical qualitative models from reaction networks.

The inputs of the tool are specified in the following files located in the folder "case_studies/model/input_files" (where 'model' denotes the name of the model, either "prozone" or "unary_vs_binary"):
- init_cond.txt: initial conditions of the model;
- reaction_rates.txt: the kinetic rates of the reactions;
- nb_sampl_interv.txt: the number of sampling intervals;
- model_reactions.txt: the reactions of the model;
- model_reactions_alphabet.txt: the alphabet of the model
- mass_invar.txt: a mass invariant of the model.

The following command:
- make 
removes any output file, compile each model and translate all dot into pdf. The generated results (both .dot and .pdf files) are stored in the folder "case_studies/model/output_files" where "model" denotes the name of the model (here either "prozone" or "unary_vs_binary").

It is also possible to decompose the computation of this script stepwise. 


The following command removes any existing output:
- make clean

The following commands generate the results from the inputs:
- make prozone: to generate the results for the prozone case study
- make unary_vs_binary: to generate the results for the unary vs binary case study
- make models: to compile all models

The following command compiles dot files into pdf files:
- make pdf











