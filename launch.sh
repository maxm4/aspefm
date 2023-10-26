# requirements: conda install -c potassco clingo=5.4.0
# conda install -c ibmdecisionoptimization cplex
# conda install -c conda-forge future
# for parse_output.py: conda install numpy; conda install pandas
# installation of mparser, an independent converter of metabolic networks to ASP and more is highly recommended

# format: (clingo clingoLP.py solveLP.lp4 [file.lp4] -c nstrict=0 -c epsilon="(1,1)" -n [EMs number] --heuristic Domain --enum-mode domRec > [output.txt])

# example
echo "== compute all EFMs of Covert Palsson 2003 model =="
clingo clingoLP.py solveLP.lp4 data/covert_palsson.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_covert_palsson.txt
echo "> written in output_covert_palsson.txt"

# efms under regulation
echo "== compute EFMs of Covert Palsson 2003 model under regulation constraints =="
clingo clingoLP.py solveLP.lp4 data/covert_palsson.lp4 data/covert_palsson_addconstr.lp4 data/covert_palsson_constr.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_covert_palsson_regul.txt
echo "> written in output_covert_palsson_regul.txt"

# biofuel
echo "== compute all EFMs of the Biofuel model  =="
clingo clingoLP.py solveLP.lp4 data/toy_model_biofuel.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_biofuel.txt
echo "> written in output_biofuel.txt"

# mcfms
echo "== compute MCFMs of Biofuel model under solution-space modifying constraints =="
clingo clingoLP.py solveLP.lp4 data/toy_model_biofuel.lp4 data/toy_model_additional_constraint.lp4 -c nstrict=0 -c epsilon="(1,1)" -n 0 -c accuracy=10 --heuristic Domain --enum-mode domRec > output_biofuel_mcfm.txt
echo "> written in output_biofuel_mcfm.txt"

# parse output
echo "== parse EFMs of Covert Palsson 2003 model into a csv file =="
python parse_output.py output_covert_palsson.txt --csv output_covert_palsson.csv
echo "> written in output_covert_palsson.csv"

