<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions that compute solutions analytically
* General Compartment Model functions that compute solutions numerically
* A flexible mechanism that handles the event schedule of clinical trials

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. There is indeed still a lot of work to do! The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC.

See the prototype manual (`torstenManual.pdf`) for more information and guidance on the examples. 

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Install
-------
Torsten is still under development and has not been added to stan-dev. The code is on github on a forked repo of stan-dev. To install cmdStan with torsten, run the shell script `script setupTorsten.sh`.

Examples
---------
For each model, we provide the following files:
* *modelName*.stan
* *modelName*.data.R
* *modelName*.init.R
* *modelName*Simulation.R 

The simulation file can be used to create the data and init files. 

A detailed description of each example can be found in the manual, but I will give a brief overview here. 

TwoCptModelExample and GenTwoCptModelExample both handle a two compartment model with first order absorption for a single patient. TwoCptModelExample computes an analytical solution, while GenTwoCptModelExample uses a numerical approximation. Given the relative simplicity of the model and the data, the fitting occurs on the order of a few seconds.  

TwoCptModelPopulation deals with the same model, with an additional level, since this time we are administering the drug to multiple patients. We model inter-individual differences in the PK parameters as random variations.

Both `TwoCptModelExampleSimulation.R` and `TwoCptModelSimulation.R` use *mrgsolve* (http://metrumrg.com/opensourcetools.html) to simulate data, and can serve as templates to simulate more PKPD data sets. I encourage interested users to experiment with parameter values and different clinical events. Changing the model is also a possibility, though this may lead to model misspecifications.  

The third example, multiDoseME2PK1, is a more elaborate model taken from the workshop *Getting Started with Bayesian PKPD Modeling Using Stan*, presented at the PAGE 2016 conference. It has the merit of illustrating how Stan and Torsten can be used to tackle more difficult problems. `MetrumBayesianStanPAGE_example3.pdf` is an excerpt from the workshop slides that describes the model. Just a heads up: Given the complexity of the model and the data, it took me about 8 hours to run this model! 

Future Development
------------------
* Add summary of output for all models in manual
* Create project directory with R scripts to facilitate simulating data, running Stan, and doing posterior predictive checks. 
