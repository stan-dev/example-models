# Bayesian Cognitive Modeling

Models from the book Bayesian Cognitive Modeling: A Practical Course (2014) by Michael Lee and Eric-Jan Wagenmakers.

Original JAGS and WinBUGS models are available here: [http://bayesmodels.com/](http://bayesmodels.com/)  
On the same address you can also find the first two parts of the book and answers to the excercises. 

Stan translation of the models was done by Martin Smira (smira.martin@gmail.com). Special thanks to Bob Carpenter and the whole Stan community.

#### NOTES 

1. All models are prepared for running from R using rstan. Model codes are included in R scripts as a string. This will change in future - models and data will be separated in standalone files that could be used in PyStan or CmdStan. We’re looking for a volunteer to pull the models out of the R code and generate data for them in a standalone file that could be used in PyStan or CmdStan.

2. We've created a google spredsheet [[here](https://docs.google.com/spreadsheets/d/1DQ_4733UxOqiP_o-HJCaQYMCV94rAvN4GQGHLCDylvM/edit?usp=sharing)] with a list of all models and with some notes concering model implementation. We weren't able to implement 3 models (out of 57 models) into Stan.

3. In the models from the first chapters, you will find notes concerning Stan related differences in model implementation. It is also advised to read chapter "Stan for Users of BUGS" in the Stan manual.

4. There is a chance that you'll find a mistake in the models, please let me know if you do. Also contact me if you have suggestions to improve the code, but keep in mind that we want to keep readability of the code high. However, please ignore the style inconsistency in the R plotting code. 

