========================
README
========================
------------
I)	THE PROGRAM
------------

This program simulates the ER+ breast cancer network in Zañudo et al. 2017 (https://www.biorxiv.org/content/early/2017/10/27/176214). For more detais on the model, see the accompanying manuscript:

Zañudo, J. G. T., & Albert, R. (2017).
A network modeling approach to elucidate drug resistance mechanisms and predict combinatorial drug treatments in breast cancer.
bioRxiv, 176214.
https://www.biorxiv.org/content/early/2017/10/27/176214

The model starts from a cancer steady state initial condition that can be primed (BCL2_T=1, BIM_T=1) or unprimed (BCL2_T=0, BIM_T=0). The model uses stochastic asynchronous updating, where a single node is updated at each update step and this node is chosen at random. The probability for a node to be selected at each time step depends only on if the node is fast (e.g. signaling events) or slow (e.g. transcriptional events), where we take fast probability taken to be 5 times faster than the slow probability (note that the real timescale difference between these type of events is much larger). The simulations are performed using 10,000 initial coniditions and simulated for 100 time steps, where each time step is defined as the average number of asychronous updates for a slow node to be updated. 

The different java classes perform different simulations:

•	NetworkSimulations - simulates the model under specified single, double, triple and no perturbations. 
•	NetworkSinglePerturbationSimulations - simulates the model in the presence of Alpelisib under all possible single node perturbations.
•	NetworkDoublePerturbationSimulations - simulates the model in the presence of Alpelisib under all possible double node perturbations.

The perturbations in each simulation are started at either time step 0 (for NetworkSinglePerturbationSimulations or NetworkDoublePerturbationSimulations) or time step 2 (NetworkSimulations).

To perform the simulations, the program uses the BooleanDynamicModeling java library (https://github.com/jgtz/BooleanDynamicModeling).

------------
II)	INSTRUCTIONS
------------

To run the program, go to the command line, navigate to the folder where the "BreastCancerModel.jar" file and the "lib" folder are located. Once there type the command:

java -jar BreastCancerModel.jar BreastCancerModel.txt

where "BreastCancerModel.txt" is the name of the TXT file with the functions of the ER+ breast cancer networ model. The "BreastCancerModel.txt" file has the following format:

"
#BOOLEAN RULES
Node1 *= Node2 or Node3
Node2 *= Node1 and Node2
Node3 *= ((not Node3 or Node4) and not Node1)
Node4 *= 1
Node5 *= 0
...
NodeN *= not Node1 or (Node1 and Node2)
"

In the above, the text before the "*=" symbol is the node name, while the text after the "*=" symbol is the Boolean function of the node.

NODE NAMES

For the node n*ames use only alphanumeric characters (A-Z,a-z), numbers (0-9) and "_". The reserved words for the program, which shouldn't be used for node names, are: "True", "False", "true", "false", "0", "1", "and", "or", and "not".

BOOLEAN FUNCTIONS

For the Boolean functions use only the node names, the logical operators "and", "or", "not", and the parentheses symbols ")" and "(". In case the Boolean function is constant, use "0" or "1", depending on the constant state of the function. The logical function does not need to be written in a disjunctive normal form; the program will take the logical form in the TXT file and transform it into its disjunctive normal form using the Quine–McCluskey algorithm.

------------
III)	OUTPUT
------------

The program will produce the following:

•	Tab separated TXT files "timecourse_X.txt" with the average trajectory of each node, where X specifies the perturbations simulated with the model (Alpelisib=1;Alpelisib=1_Everolimus=1;Alpelisib=1_MCL1=0;Alpelisib=1_Palbociclib=1_Fulvestrant=1;Alpelisib=1_PIM=1).
•	A tab separated TXT file "BreastCancerSinglePerturbations.txt" with the results of the simulating the model under each possible single node perturbation. The file contains the average Apoptosis (Apoptosis, Apoptosis_2, Apoptosis_3, Apoptosis_norm) and Proliferation (Proliferation, Proliferation_2, Proliferation_3, Proliferation_4, Proliferation_norm) values at the end of the simulation of each perturbation.
•	A tab separated TXT file "BreastCancerDoublePerturbations.txt" with the results of the simulating the model under each possible double node perturbation. The file contains the average Apoptosis (Apoptosis, Apoptosis_2, Apoptosis_3, Apoptosis_norm) and Proliferation (Proliferation, Proliferation_2, Proliferation_3, Proliferation_4, Proliferation_norm) values at the end of the simulation of each perturbation.

------------
IV)	SOFTWARE USED AND LICENSES
------------

JGraphT

Several functions from the JGraphT java class library by Barak Naveh and Contributors are used (https://github.com/jgrapht/jgrapht). JGraphT is available under GNU LESSER GENERAL PUBLIC LICENSE Version 2.1.

Quine-McCluskey_algorithm

An implementation of the Quine-McCluskey_algorithm in the “Term.java” and “Formula.java” classes were retrieved in 2013 from http://en.literateprograms.org/Quine-McCluskey_algorithm_(Java)?action=history&offset=20110925122251. The “Term.java” and “Formula.java” classes are available under the MIT License.

------------
V)	COPYRIGHT
------------

The MIT License (MIT)

Copyright (c) 2017 Jorge G. T. Zañudo and Réka Albert.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
