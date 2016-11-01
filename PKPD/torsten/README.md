<b> Torsten </b> is a library of C++ functions that support applications of Stan in Pharmacometrics. The current prototype provides:
* One and Two Compartment Model functions that compute solutions analytically
* General Linear Compartment Model function that computes a matrix exponential solution
* General Compartment Model functions that compute solutions numerically
  * Runge-Kutta 4th/5th order for non-stiff ODE systems
  * Backward Differentiation for stiff ODE systems
  
** This prototype is still under development ** and has been uploaded to facilitate working with the community of Stan developers. The current version was written by Charles Margossian, Bill Gillespie, and Metrum Research Group, LLC. We have recieved extensive help from the Stan development team.

See the user manual (`torstenManual.pdf`) for more information and guidance on the examples. If you have any questions, please raise an issue on GitHub or send me an e-mail at charlesm@metrumrg.com. 

Licensing
---------
The Torsten library is licensed under the BSD 3-clause license. 


Install
-------
To install cmdStan with torsten, run the shell script script `setupTorsten.sh`.

We are working with Stan's development team to create a system to add and share Stan packages. In the mean time, users can download a forked version of Stan with Torsten from GitHub. Torsten uses features from stan-dev, the development version of Stan. This means it currently requires the CmdStan interface which is the only interface up to date with stan-dev.

When Stan releases version 2.13 (expected end of November 2016), I will update Torsten to be compatible with the new realease and usable with all its interfaces, including RStan.


Examples
---------
For each model, we provide the following files:
* *modelName*.stan
* *modelName*.data.R
* *modelName*.init.R
* *modelName*Simulation.R 

The simulation file can be used to create the data and init files. 

A description of these examples and instructions on how to run them can be found in the manual. 

Under the R directory, we provide tools to run the examples via cmdStan and look at diagnostic plots.

C++ Code
--------
This C++ code for Torsten can be found on the following repos:

Math Library: https://github.com/charlesm93/math/tree/feature/issue-314-torsten-pharmacometrics

Stan Grammar: https://github.com/charlesm93/stan/tree/feature/issue-1953-exposing-torsten-pharmacometrics

Updates
-------
10/31
* Add Linear Compartment model function.
* Revise the user's manual.
* Torsten is now compatible with a development version, post v 2.12. 

08/02 
* Update Stan and Math branch to match stan-dev. This is important because of the recent bug fix in stan 2.10. 

