# dyn_qual_compil

Author: Wassim Abou-Jaoudé
Date:30 April 2015
------------------------------------------------------------------------------------------

Matlab program for the automatic generation of dynamical qualitative models from reaction networks.

The inputs of the tool are specified in the following files located in the folder "\tool\input\model" (where 'model' denotes the name of the model, either 'prozone' or 'rescaling'):
- init_cond.txt: initial conditions of the model;
- reaction_rates.txt: the kinetic rates of the reactions;
- nb_sampl_interv.txt: the number of sampling intervals;
- model_reactions.txt: the reactions of the model;
- model_reactions_alphabet.txt: the alphabet of the model
- mass_invar.txt: a mass invariant of the model.

To generate the results from the inputs, launch in Matlab the command line:

	frontend('model')

where 'model' denote the name of the model (either 'prozone' or 'rescaling').

The generated output is a .dot file stored in the folder "\tool\result".

Then one can use a visualisation software (Graphviz for example) to generate the graph from the dot file.







