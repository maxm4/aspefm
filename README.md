# aspefm

Computing Elementary Flux Modes (EFM) with Answer Set Programming (ASP)

## How to use :

This tool uses `clingo[LP]` (clingoLP.py) and the rules specified in `solveLP.lp4` to compute the EFMs of an input metabolic network.

The metabolic network should be in ASP format.

Example executions of this tool are found in the bash file `./launch.sh`.

## Required :

### Python version :

This module requires python version 3.

### All uses :

clingo version 5.4.0, compilation of clingo with Python is mandatory

	conda install -c potassco clingo

### For compatibility with Python 2 purposes :

	conda install -c conda-forge future
	
### Python installs :

pandas, numpy, scipy, sympy modules

	conda install pandas
	conda install numpy
	conda install scipy
	conda install sympy

### Clingo[LP] :

Clingo[LP] is a clingo linear programming extension, compatible with the cplex solver.

If you don't have access to a cplex license, you may install the anaconda distribution :

	conda install -c ibmdecisionoptimization cplex

Then, this extension only requires to use the script `clingoLP.py`, which is provided on this repository.

Here is how to use aspefm from a metabolic network file encoded in ASP :

	clingo clingoLP.py [FILE] solve[LP].lp4 -c nstrict=0 -n 0 -c accuracy=10

### Minimization heuristics :

Minimization heuristics are working with clingo extensions `clingo[LP]`:

To use them, add the following options when running `clingo` :

	--heuristic Domain --enum-mode domRec


### Constraints :

You can constraint the metabolic network by specifying for example that you only want elementary modes containing the biomass reaction.

ASP constraints should be expressed this way:

	:- not support("biomass"). 
	
This imposes the solver to have the biomass reaction in the solution.

Alternatively,

	:- support("biomass"). 
	
Imposes the solver to not compute the biomass reaction in the solution.
