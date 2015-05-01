# dyn_qual_compil

Author: Wassim Abou-Jaoudé
Date:30 April 2015
------------------------------------------------------------------------------------------


We need few days to organize the repository so that the scripts work properly. Sorry about it. 


The following instructions will work soon.

Matlab program for the automatic generation of dynamical qualitative models from reaction networks.

The inputs of the tool are specified in the following files located in the folder "case_studies/model/input_files" (where 'model' denotes the name of the model, either 'prozone' or 'unary_vs_binary'):
- init_cond.txt: initial conditions of the model;
- reaction_rates.txt: the kinetic rates of the reactions;
- nb_sampl_interv.txt: the number of sampling intervals;
- model_reactions.txt: the reactions of the model;
- model_reactions_alphabet.txt: the alphabet of the model
- mass_invar.txt: a mass invariant of the model.

To generate the results from the inputs, use the command line:
- make prozone: to generate the results for the prozone case study
- make unary_vs_binary: to generate the results for the unary vs binary case study
- make models: to compile all models

The generated results are .dot files which are stored in the folder "case_studies/model/output_files" where "model" denotes the name of the corresponding case study (here either "prozone" or "unary_vs_binary").

Then use the command line: make pdf, to generate the resulting graphs in .pdf format from the .dot file.







