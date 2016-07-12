<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions, that compute solutions analytically
* General Compartment Model functions, that compute solutions numerically
* A flexible mechanism that handles the event schedule of clinical trials

This prototype is still under development and has been uploaded to facilitate working with the community of Stan developers. The current prototype was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC.

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 

Examples
---------
These examples are meant to illustrate ways in which Torsten can be used. This directory contains the Stan model and the data files for each example. In addition the project directory provides a working structure to simulate PKPD data (using the R package mrgsolve) and test it against Stan models. Torsten is still under development and does not yet appear in stan-Dev, but the code is available on GitHub (on a forked version of stan-dev). The bash file, setupTorsten.sh, downloads and installs cmdstan-dev and the Torsten library.  
