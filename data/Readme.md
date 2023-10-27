# Example models

Please use the module *[mparser](https://github.com/maxm4/mparser)* for creation of metabolic networks for usage with *aspefm*.

```covert_palsson.lp4```:

```covert_palsson_constr.lp4```:

```covert_palsson_addconstr.lp4```:

Examples and regulation constraints from :

- Covert and Palsson, 2003, "Constraints-based models: Regulation of Gene Expression Reduces the Steady-state Solution Space".

```toy_model_biofuel.lp4```:

- example model for usage with F-A-M-E (maintained by Brett Olivier),

```toy_model_additional_constraint.lp4```:

- specifies that ```fuel1```, ```fuel2```, ```biomass``` are positive inputs at the same time,
- generates solutions that are not EFMs if not filtered out.

```toy_model_constraint.opt```:

- table of table of tuples describing the two following constraints:

- ```-R09 -R15 -R17 <= 0```

- ```R02 -R03 <= 0```

```toy_model_biofuel.efmc```:

Pickled data structure containing the stocihiometric matrix of the previous network for use with extension *EFMChecker*.

Code for constructing this structure distributed with *mparser*, simplifications will be done at a later date.

```mcs_example_hadicke.lp4```:

```mcs_example_target.lp4```:

```mcs_example_constraint_1.lp4```:

```mcs_example_constraint_2.lp4```:

```mcs_example_wanted_constraint.lp4```:

Dual network for MCSs computation with formulation from:

- von Kamp and Klamt, 2014, "Enumeration of smallest intervention strategies in genome-scale metabolic networks".

Network and constraints adapted from:

- Hadicke and Klamt, 2011, "Computing complex metabolic intervention strategies using constrained minimal cut sets".

```mcs_example.xml```:

SBML Model of the previous network for use with extension *MCSChecker*.
