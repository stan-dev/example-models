# setupTorsten
# v 1.1
# Run this bash file in the directory where you want to install cmdstan
# with the Torsten functions. 
#
# update 1.1: download last version of cmdStan/dev with which Torsten
# was tested.  

#!/bin/bash
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
git reset --hard 53b4041
git clone https://github.com/charlesm93/stan.git
cd stan
git checkout feature/issue-1953-exposing-torsten-pharmacometrics
cd lib
rm -r stan_math
git clone https://github.com/charlesm93/math.git
mv math stan_math
cd stan_math
git checkout feature/issue-314-torsten-pharmacometrics
