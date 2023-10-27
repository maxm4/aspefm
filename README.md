# aspefm

Computing Elementary Flux Modes (EFM) with Answer Set Programming (ASP).

## Context :

Elementary Flux Modes are minimal sets of enzymes that operate at steady state with all irreversible reactions proceeding in the appropriate direction.

The enumeration of EFMs is a difficult task. It requires the resolution of combinatorial problems on metabolic networks, and the integration of appropriate biological constraints to help calculations.

We propose to use the SAT-based power of ASP constraint logic programming resolution to reduce the hurdle of obtaining pathways of interest with EFMs on large-scale networks.

## How to use :

This tool uses `clingo[LP]` (*clingoLP.py*) and the Answer Set Programming (ASP) rules specified in `solveLP.lp4` to compute the Elementary Flux Modes (EFMs) of a given metabolic network.

The input metabolic network should be in ASP format.

The module *[mparser](https://github.com/maxm4/mparser)* should be used for converting metabolic networks into ASP rules.

Example executions of the tool are found in the bash file `./launch.sh`, `./ext_launch.sh`, `./mcs_launch.sh`.

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

### clingo[LP] :

[clingo[LP]](https://github.com/potassco/clingoLP/) is a linear programming extension for the solver [clingo](https://github.com/potassco/clingo/), compatible with the [IBM cplex](https://www.ibm.com/products/ilog-cplex-optimization-studio) solver.

If you don't have access to a cplex license, you may install the anaconda distribution :

	conda install -c ibmdecisionoptimization cplex

Then, this extension only requires to use the script `clingoLP.py`, which is provided on this repository.

Here is how to use *aspefm* from a metabolic network file encoded in ASP :

	clingo clingoLP.py [FILE] solve[LP].lp4 -c nstrict=0 -n 0 -c accuracy=10

### Minimization heuristics :

Minimization heuristics are working with clingo extensions `clingo[LP]`:

To use them, add the following options when running `clingo` :

	--heuristic Domain --enum-mode domRec


### Constraints :

One can constrain the metabolic network by specifying for example that only the elementary modes containing the biomass reaction are needed.

ASP constraints should be expressed this way:

	:- not support("biomass"). 
	
This imposes the solver to have the biomass reaction in the solution.

Alternatively,

	:- support("biomass"). 
	
Imposes the solver to not compute the biomass reaction in the solution.

### Acknowledgements :

If you want to use *aspefm*, we refer you to our paper: [Mahout, M., Carlson, R. P. and Peres, S. Answer Set Programming for Computing Constraints-Based Elementary Flux Modes: Application to Escherichia coli Core Metabolism. Processes 8, 1649 (2020)](https://doi.org/10.3390/pr8121649).

Additionally, `clingo[LP]` is property of the University of Potsdam, Germany and collaborators.

Please refer to the following publication: [Janhunen, T. et al. "Clingo goes linear constraints over reals and integers." CoRR abs/1707.04053, (2017)](https://doi.org/10.1017/S1471068417000242)).


