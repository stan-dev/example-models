<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions, that compute solutions analytically
* General Compartment Model functions, that compute solutions numerically
* A flexible mechanism that handles the event schedule of clinical trials

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC.

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Install
-------
Torsten is still under development and has not been added to stan-dev. The code is on github on a forked repo of stan-dev. To install cmdStan with torsten, run the shell script setupTorsten.sh.

Examples
---------
These examples illustrate how Torsten can be used. The two models, TwoCptModelExample and GenTwoCptModelexample, deal with a Two Compartment model with a first order absorption. The first example uses an analytical solution, the latter a numerical one. TwoCptModelSimulation.R simulates PKPD data using the R package mrgsolve (http://metrumrg.com/opensourcetools.html). I encourage users to create different data sets with this template, and fit them. 

More examples to come, along with summaries of output and documentation.  
