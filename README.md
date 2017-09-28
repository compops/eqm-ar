eqm-ar
======

Equalisation Maximisation (EqM) and Expectation Maximisation (EM) for AR(p) processes

This code was downloaded from < https://github.com/compops/eqm-ar/ > and contains the code used to produce the results in the poster presentation

J. Dahlin, F. Lindsten and T. B. Schön, **Parameter inference in AR processes with missing data**. ERNSI Workshop, Maastricht, NL, 2012.

which is available as from < http://arxiv.org/pdf/1308.4601 >.

Included files
--------------

**RUNME**
The code reproduces the comparision presented in the right column in the result section of the poster.

Helpfiles and subroutines are found in the folder *helpers*:
- ARemsub:		EM algorithm for AR(p) process with missing data
- AReqmsub:		EqM algorithm for AR(p) process with missing data
- ARstdsub:		LS for AR(p) process with missing data
- runARmodel:		Generates data from random AR(p) systems with missing observations
- BuildPhi:		Constructs the Hankel matrix of past inputs and outputs.
