# Multithreading and Map-Reduce in Stan 2.18.0: A Minimal Example

v1 - Richard McElreath - 19 Sep 2018

## Introduction

In version 2.18.0, Stan can now use parallel subthreads to speed computation of log-probabilities. So while the Markov chain samples are still necessarily serial, running on a single core, the math library calculations needed to compute the samples can be parallelized across different cores. If the data can be split into pieces ("shards"), then this can produce tremendous speed improvements on the right hardware.

The benefits are not automatic however. You'll have to modify your Stan models first. You might also have to do a little bit of compiler configuration. Finally, the 2.18.0 is only available for ``cmdstan``, the command line version of Stan, at this moment. It should be available for ``rstan`` and other interfaces soon. So if you haven't used cmdstan before, you'll need to adjust your workflow a little.

Here I present a minimal example of modifying a Stan model to take advantage of multithreading. I also show the basic workflow for using ``cmdstan``, in case you are used to ``rstan``.

## Overview

The general approach to modifying a Stan model for multithreading has three steps.

(1) You must repack the data and parameters into slices, known as "shards". This looks weird the first time you do it, but it really is just repacking. Nothing fancy is going on.

(2) You must rewrite the model to use a special function to compute the log-probability for each slice ("shard") of the data. Then in model block, you use a function called `map_rect` to call this function and pass it all the shards. 

(3) Configure your compiler to enable multithreading. Then you can compile and sample as usual.

In the rest of this tutorial, I'll build up a simple first without multithreading. This will serve to clarify the changes and show you how to work with `cmdstan`, if you haven't before. Then I'll modify the data and model to get multithreading to work.

If you use Windows, it looks like multithreading isn't working yet (see <https://github.com/stan-dev/math/wiki/Threading-Support>). So flip over to a Linux partition, if you have one. If you don't, hope is that an upcoming update to RTools will make things work.

## Prepping the data

Let's consider a simple logistic regression, just so the model doesn't get in the way. We'll use a reasonably big data table, the football red card data set from the recent crowdsourced data analysis project (<https://psyarxiv.com/qkwst/>). This data table contains 146,028 player-referee dyads. For each dyad, the table records the total number of red cards the referee assigned to the player over the observed number of games. 

The ``RedcardData.csv`` file is provided in the repository here. Load the data in R, and let's take a look at the distribution of red cards:
```R
d <- read.csv( "RedcardData.csv" , stringsAsFactors=FALSE )
table( d$redCards )
```
```text
     0      1      2 
144219   1784     25
```
The vast majority of dyads have zero red cards. Only 25 dyads show 2 red cards. These counts are our inference target.

The motivating hypothesis behind these data is that referees are biased against darker skinned players. So we're going to try to predict these counts using the skin color ratings of each player. Not all players actually received skin color ratings in these data, so let's reduce down to dyads with ratings:
```R
d2 <- d[ !is.na(d$rater1) , ]
out_data <- list( n_redcards=d2$redCards , n_games=d2$games , rating=d2$rater1 )
out_data$N <- nrow(d2)
```
This leaves us with 124,621 dyads to predict. 

At this point, you are thinking: "But there are repeat observations on players and referees! You need some cluster variables in there in order to build a proper multilevel model!" You are right. But let's start simple. Keep your partial pooling on idle for the moment.

Now to use these data with ``cmdstan``, we need to export the table to a file that ``cmdstan`` can read in. This is made easy:
```R
library(rstan)
stan_rdump(ls(out_data), "redcard_input.R", envir = list2env(out_data))
```

## Making the model

A Stan model for this problem is just a simple logistic (binomial) GLM. I'll assume you know Stan well enough already that I can just plop the code down here. It's contained in the ``logistic0.stan`` file in the repository. It's not a complex model:
```c++
data {
  int N;
  int n_redcards[N];
  int n_games[N];
  real rating[N];
}
parameters {
  vector[2] beta;
}
model {
  beta ~ normal(0,1);
  n_redcards ~ binomial_logit( n_games , beta[1] + beta[2] * to_vector(rating) );
}
```
If you were going to run this model in ``rstan``, you'd enter now:
```R
m0 <- stan( file="logistic0.stan" , data=out_data , chains=1 )
```
But we're going to instead build and run this in ``cmdstan``.

## Installing and Building cmdstan

You'll need to download the 2.18.0 `cmdstan` tarball from here: <https://github.com/stan-dev/cmdstan/releases/tag/v2.18.0>. Also grab the 2.18.0 user's guide while you are at it. It has a nice new section on this topic.

If you are used to unix and to `make` workflows, you can expand the thing and build it now. But if you are not, don't worry. I'm here to help. Many of my colleagues are very capable scientists who have never done much work in a unix shell. So I know even highly technical people who need some basic help here. I am going to assume, however, that you aren't using Windows. Because if you are, you can install and use `cmdstan` fine, but the multithreading won't work (yet).

Begin by putting the `cmdstan-2.18.0.tar.gz` file wherever you want `cmdstan` to live. You can move it later. Then open a terminal and navigate to that folder. Decompress with:
```bash
tar -xzf cmdstan-2.18.0.tar.gz
```
Enter the new `cmdstan-2.18.0` folder and build with:
```bash
make build
```
Your compiler will toil for a bit and then signal its satisfaction.

## Building and Running the Basic Model in cmdstan

Make a new folder called `redcard` inside your `cmdstan-2.18.0` folder and drop the `redcard_input.R` and `logistic0.stan` files in it. To build the Stan model, while still in the `cmdstan-2.18.0` folder, enter:
```bash
make redcards/logistic0
```
Again, compiler toil. Satsifaction. To sample from the model:
```bash
cd redcard
time ./logistic0 sample data file=redcard_input.R
```
I added the `time` command to the front so that we get a processor time report at the end. We'll want to compare this time to the time you get later, after you enable multithreading. On my machine, I get:
```
real    2m22.694s
user    2m22.369s
sys     0m0.247s
```
The first line is the "real" time, the one you probably care about for benchmarking.

By default, the output file is called `output.csv`. You'll want to go back into R to analyze it. Just open R in the same folder and enter:
```R
library(rstan)
m0 <- read_stan_csv("output.csv")
```
Now you can proceed as usual to work with the samples.

## Rewriting the Model to Enable Multithreading

The first thing to do in order to get multithreading to work is to rewrite model. Our rewrite needs to accomplish two things. 

First, it needs to restructure the data into multiple "shards" of the same length. Each shard can be sent off into its own thread. The right number of shards depends upon your data, model, and hardware. You could make each dyad, in our example, into a shard. That would be pretty inefficient however, because vectorizing the calculation is usually beneficial, and you probably don't have enough processor cores to handle that many parallel threads at once. At the other end, two shards is better than one, and three is nearly always better than two. At some point, we reach an optimal number that trades off efficiencies within each shard with the ability to have parallel calculations.

Second, it needs to use Stan's new `map_rect` function to update the log-probability. This will mean moving the log-probability calculation to a special function and then calling that function in the model block using `map_rect`.

We'll take on each of these in turn.

### Making Shards

In this example, I am going to make 7 shards of equal size. I didn't find this number by any reasoning. It is just that we have 124,621 dyads, and that divides nicely by 7 into 17,803 dyads per shard. So 7 shards of 17,803 dyads each.

You'll want to do the splitting into shards in the `transformed data` block of your Stan model. I've built a modified version of `logistic0.stan` and named it, creatively, `logistic1.stan`. It's in the repository. Here is its `transformed data` block:
```c++
transformed data {
  // 7 shards
  // M = N/7 = 124621/7 = 17803
  int n_shards = 7;
  int M = N/n_shards;
  int xi[n_shards, 2*M];  // 2M because two variables, and they get stacked in array
  real xr[n_shards, M];
  // an empty set of per-shard parameters
  vector[0] theta[n_shards];
  // split into shards
  for ( i in 1:n_shards ) {
    int j = 1 + (i-1)*M;
    int k = i*M;
    xi[i,1:M] = n_redcards[ j:k ];
    xi[i,(M+1):(2*M)] = n_games[ j:k ];
    xr[i] = rating[ j:k ];
  }
}
```
The way that shards work is that we have to pack all of the variables particular to each shard into three arrays: an array of integer data, an array of real (continuous) data, and an array of parameters. I'm going to call the array of integers `xi`, the array of reals `xr`, and the array of parameters `theta`. 

First, consider the array of reals, `xr`. It's the simplest. We need to pack into this the `rating` values. So we define:
```
real xr[n_shards, M];
```
The value `M` is just the number of dyads per shard. Then in the loop at the bottom of the block, we pack in the `rating` values from each shard, iterating over the indexes.

The array of integers is slightly harder. We have two integer variables to pack in: `n_redcards` and `n_games`. But these need to be stacked together in a single vector. The code above defines therefore:
```
int xi[n_shards, 2*M];
```
Each row of `xi` is twice as long as `xr`, because we need two variables. These get packed in the loop at the bottom. We will unpack them later.

The array of parameters works the same way. But in this example, it is super easy. There are no parameters special to each shard, because the parameters are global (no varying effects, for example). So we can just define a zero-length array for each shard:
```
vector[0] theta[n_shards];
```
If there were actually parameters to pack in here, we'd define them in either the `parameters` or `transformed parameters` block.

Hopefully that explains the shard construction. Your data will need their own packing, but the principles are the same. 

### Using map_rect

In the new model, we replace the model block with:
```c++
model {
  beta ~ normal(0,1);
  target += sum( map_rect( lp_reduce , beta , theta , xr , xi ) );
}
```
Gone is the call to `binomial_logit`. Instead we call `map_rect` and give it five arguments:

(1) `lp_reduce`: This is name of our new function (we'll write it soon) that computes the log-probability of the data.

(2) `beta`: This is a vector of global parameters that do not vary among shards. In this case, it is just the vector defined in the `parameters` block as before. Nothing changes.

(3) `theta`: This is the array of parameters that vary among shards. We defined this in the previous section.

(4) `xr` and (5) `xi`: These are the data arrays we defined above.

Then we make a new block, at the very top of our Stan model, to hold the function `lp_reduce`:
```c++
functions {
  vector lp_reduce( vector beta , vector theta , real[] xr , int[] xi ) {
    int n = size(xr);
    int y[n] = xi[1:n];
    int m[n] = xi[(n+1):(2*n)];
    real lp = binomial_logit_lpmf( y | m , beta[1] + to_vector(xr) * beta[2] );
    return [lp]';
  }
} 
```
You can name this function whatever you want. But it has the have arguments of types `vector`, `vector`, `real`, and `int`. Those match the arguments in the call to `map_rect`.

Inside the function, we first unpack the data into separate arrays again. Then the familiar call to `binomial_logit`, but this time storing the result in `lp`. Then we convert `lp` to a vector and return it. 

When you run the model, this vector lands back in the model block, `map_rect` collects all the different returns from each shard, sums them together, and updates the log-probability. Samples are born.

## Running the Multithreaded Model

### Configure Compiler

Now we're ready to go. Almost. First you need to add some compiler variables. Go to your `cmdstan-2.18.0` folder from earlier. Then enter in the shell:
```bash
echo "CXXFLAGS += -DSTAN_THREADS" > make/local
```
What this does is add a special flag for `cmdstan` when it compiles your Stan models. If you are using g++ (instead of clang), you will also likely need one more flag:
```bash
echo "CXXFLAGS += -pthread" >> make/local
```
On my Linux machine, that was enough to make it work. On my Macintosh, all that was needed was `-DSTAN_THREADS`. If you want to see the `local` file you just made, do:
```bash
cat make/local
```

### Build and Run

Now build the Stan model:
```bash
make redcard/logistic1
```
And finally we can switch into the `redcard` folder, set the number of threads we want to use (-1 means all), and go:
```bash
cd redcard
export STAN_NUM_THREADS=-1
time ./logistic1 sample data file=redcard_input.R
```
The timings that I get for this model are:
```
real  1m42.789s
user  5m10.812s
sys   0m44.104s
```
So that's 1m43s versus 2m23s from before. We saved less than a minute by multithreading. Not that big a deal, in absolute time. But as a proportion it is really something. And once you start adding varying effects to this model, or have a larger data table, that's quite a speedup you can expect.

## More

There is a new section in the 2.18 user manual that contains examples of using `map_rect`. Get it here: <https://github.com/stan-dev/cmdstan/releases/tag/v2.18.0> and see page 237. It also contains an example of a hierarchical model.
