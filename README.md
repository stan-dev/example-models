## Example Models

This repository holds open source Stan models, data simulators, and real data.  There are models translating those found in books, most of the BUGS examples, and some basic examples used in the manual.

#### Books

* [Applied Regression Modeling](https://github.com/stan-dev/example-models/wiki/ARM-Models) (Gelman and Hill 2007)
* [Bayesian Cognitive Modeling](https://github.com/stan-dev/example-models/tree/master/Bayesian_Cognitive_Modeling) (Lee and Wagenmakers 2014)
* [Bayesian Population Analysis using WinBUGS: A Hierarchical Perspective](https://github.com/stan-dev/example-models/tree/master/BPA) (Kéry and Schaub 2012)

#### BUGS examples

- [BUGS Models (volumes 1--3)](https://github.com/stan-dev/example-models/wiki/BUGS-Examples)

#### Basics

- Basic Distributions
- Basic Estimators

## Organization

Each example model should be organized to include the following

1.  Stan program(s) implementing the model (and variants),
2.  program to simulate data (in R or Python),
3.  simulated data file itself (for now in .data.R dump format). For each Stan file(s), such as foo.stan and bar.stan, there must be a foo.data.R file and a bar.data.R file in the same subdirectory as the .stan file (even if foo.data.R is the same as bar.data.R).
4.  summary of output fit for simulated data file (text),
    (a) check diagnostics for post-warmup iterations: n_divergent is zero for all post-warmup iterations, tree_depth < max_depth for all post-warmup iterations,
    (b) check for all parameters (and unnormalized density lp__) that R_hat < 1.1 and n_eff > 20, and
    (c) if diagnostics fail improve model with a note explaining the fix.
5.  any real data sets with summary of fit (in original and Stan-readable format, with munging program from original to Stan-readable),
    (a) check diagnostics (as above),
    (b) check convergence (as above), and
    (c) if fit is poor then improve model with a note explaining the fix.
6.  citations to source for model and/or data,
7.  documentation needed to understand the model (LaTeX, text, or HTML), including discussion of parameterizations and/or optimizations and their relations to performance (and to each other if there are multiple models),
8.  keywords or other tags to help organize by category (e.g., from manual, from BUGS volume, from book, involving logistic regression, application area such as population model, model type such as IRT or mark-recapture, etc.)
9.  author/copyright-holder info and open-source license info if not new BSD.

The idea is to encourage people to go through *all* of these steps for their models, particularly 3 and 4, which often get overlooked. And if the example cannot be executed via rstan::stan_demo(), then something is wrong.

### Licensing

All of the example models are copyrighted by their author(s) or assignees under the new BSD license unless another open-source license is explicitly stipulated in the directory containing the model.
