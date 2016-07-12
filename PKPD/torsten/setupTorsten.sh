# setupTorsten
# v 1.0
# Run this bash file in the directory where you want to install cmdstan
# with the Torsten functions. 

#!/bin/bash
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
git clone https://github.com/charlesm93/stan.git
cd stan
git checkout feature/issue-1953-exposing-torsten-pharmacometrics
cd lib
rm -r stan_math
git clone https://github.com/charlesm93/math.git
cd math
git checkout feature/issue-314-torsten-pharmacometrics
cd ..
mv math stan_math