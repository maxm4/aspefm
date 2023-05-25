# requirements: conda install -c potassco clingo=5.4.0
# conda install -c ibmdecisionoptimization cplex
# conda install -c conda-forge future
# for parse_output.py: conda install numpy; conda install pandas

# format: (clingo clingoLP.py solveLP.lp4 [file.lp4] -c nstrict=0 -c epsilon="(1,1)" -n [EMs number] --heuristic Domain --enum-mode domRec > [output.txt])

# example
clingo clingoLP.py solveLP.lp4 data/covert_palsson.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_covert_palsson.txt

# efms under regulation
clingo clingoLP.py solveLP.lp4 data/covert_palsson.lp4 data/covert_palsson_addconstr.lp4 data/covert_palsson_constr.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_covert_palsson_regul.txt

# biofuel
clingo clingoLP.py solveLP.lp4 data/toy_model_biofuel.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_biofuel.txt

# mcfms
clingo clingoLP.py solveLP.lp4 data/toy_model_biofuel.lp4 data/toy_model_additional_constraint.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_biofuel_mcfm.txt

# parse output
python parse_output.py output_covert_palsson.txt --csv output_covert_palsson.csv


