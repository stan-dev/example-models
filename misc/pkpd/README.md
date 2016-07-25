# Draft of PKPD library

WARNING: THESE FILES ARE IN A DRAFT STATE!

This means: Code has been tested and is expected to work
correctly. However, tests have not yet been formalized nor automated
for the involved functions. Moreover, the performance of the codes may
still be sub-optimal and are subject to further optimization.

## Introduction

These files are an attempt to ease writing PKPD models in
Stan. Pharmacokinetic (PK) models describe the relationship of drug
concentration in a given subject over time as a function of the dosing
history. PK models facilitate mass action kinetics laws to describe
the drug absorption, distribution and elimination process with a
compartmental approximation.

## Dosing History

The dosing events are discontinuous changes of system states. For
example, an oral route of administration translates into a first
absorption from a dosing compartment into the main compartment. At
each dosing event the specified dose is added to the dosing
compartment such that the amount changes discontinuously.

As such the dosing history is an important part of the data and must
be handled accordingly. As PK systems can be described by ordinary
differential equations (ODE), this implies important properties of the
solution. In particular, defining a reference state y0 at time t0
fully defines the solution at finite times, given the solution
operator U(t0,t) of the ODE system. As reference state it is
convenient to choose the last dosing event which has occurred and
forward develop the solution to finite times using the solution
operator. The solution operator can be given as analytic solution or
as result from an ODE integrator.

At each event the ODE system is stopped, the event applied and then
the solution is continued to be calculated.

## Regular Dosing History

However, in the treatment of patients we often encounter a regular
dosing pattern such as dose $d$ every $\tau$ time units (24h for daily
dosing). Luckily, these dosing patterns can be efficiently handled for
analytically solvable systems and hence represent an important special
case. Analytically solvable systems are of the form of a matrix
exponential exp(At) and the linear super-position principle holds. As
consequence a regular dosing situation can quickly solved using the
geometric series.

For ODE based PK systems the above simplifications cannot be applied
and as such regular dosings patterns are not supported by the library
for ODE problems.

## NONMEM Data Sets

A widely known data format to describe the longitudinal data of each
patient is the NONMEM format. As such convenient functions are
provided to map a NONMEM data set to the internal format.

## Usage

An example program is provided with `oral_1cmt.stan`. The program
starts with including the utility functions. The user only has to
define the functions `pk_system` and `pk_system_addl` which receive
the parameters being sampled and must feed these into the PK model of
choice which is effectively the solution operator U(t0,t) for the
model at hand.

Currently, the library includes a single analytical model in two
variants and an example ODE based Michaelis Menten model is
provided. The full analytical model has the structure of three
compartments. Compartment 1 transfers mass with rate k12 to cmt 2 and
mass is eliminated from cmt 1 with rate k1. Compartment 2 receives
mass with rate k12 from cmt 1 and mass exits with rate k2 which flows
into compartment 3. The model with all three compartments is the
`pk_1cmt_metabolite_depot` and an additional variant lacks the depot
cmt. This allows to cover the situations

- 1cmt oral dosing with calculation of the AUC (set k1=k12)
- IV dosing and metabolite with AUC calc
- In fact, with appropriate super-position an arbitrary amount of
  metabolites can be computed from this; just use the given function
  multiple times.

Since the include functionality is used which is supported by rstan
stanc_builder function, an intermediate Stan file `oral_1cmt_run.stan`
is created.

# Files

- `oral_1cmt.stan` example Stan program (1cmt oral)
- `oral_1cmt.R` runs the example Stan program with simulated data
- `oral_1cmt_sim.R` simulate (somewhat realistic) data scenario
- `oral_1cmt_run.data.R` simulated data in Stan format
- `oral_1cmt_run_nm.csv` simulated data in NONMEM format
- `oral_1cmt_run.stan` Stan program with includes (DO NOT CHANGE, gets
overwritten)

- `oral_1cmt_mm.stan` example Stan program (1cmt oral, MM elimination, ODE)
- `oral_1cmt_mm.R` runs the example Stan program with simulated data
- `oral_1cmt_mm_sim.R` simulate (somewhat realistic) data scenario
- `oral_1cmt_mm_run.data.R` simulated data in Stan format
- `oral_1cmt_mm_run_nm.csv` simulated data in NONMEM format
- `oral_1cmt_mm_run.stan` Stan program with includes (DO NOT CHANGE, gets
overwritten)

- `model_lib.stan` model library, Stan user functions
- `utils.stan` various utility Stan user functions
