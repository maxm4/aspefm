clingoLP="clingo clingoLP.py"
clingoLPParams=(-c nstrict=0 -c epsilon="(1,2)" -c accuracy=10)
clingoLPParams="${clingoLPParams[@]}"
minSupport="--heuristic Domain --enum-mode domRec"
clingoParams="--stats=2"

echo "== only compute EFM solutions when with positive constraints using EFMChecker  =="
$clingoLP $clingoLPParams $clingoParams $minSupport solveLP.lp4 extensions/efm_checker.py -c efmcheckfile=\"data/toy_model_biofuel.efmc\" data/toy_model_biofuel.lp4 data/toy_model_additional_constraint.lp4 -n 0 > output_biofuel_efms.txt
echo "> written in output_biofuel_efms.txt"

echo "== debug EFM computation  =="
$clingoLP $clingoLPParams $clingoParams $minSupport solveLP.lp4 extensions/debug_propagator.py -c debug=2 -c trace=1 extensions/efm_checker.py -c efmcheckfile=\"data/toy_model_biofuel.efmc\" data/toy_model_biofuel.lp4 data/toy_model_additional_constraint.lp4 -n 0 > output_biofuel_efms_debug.txt 2> output_biofuel_efms_trace.txt
echo "> written in output_biofuel_efms_debug.txt and output_biofuel_efms_trace.txt"

echo "== profile EFM computation  =="
$clingoLP $clingoLPParams $clingoParams $minSupport solveLP.lp4 extensions/profiler_extension.py extensions/efm_checker.py -c efmcheckfile=\"data/toy_model_biofuel.efmc\" data/toy_model_biofuel.lp4 data/toy_model_additional_constraint.lp4 -n 0 > output_biofuel_efms_profile.txt
echo "> written in output_biofuel_efms_profile.txt"

echo "== remove solutions thanks to second call to LP  =="
$clingoLP $clingoLPParams $clingoParams $minSupport solveLP.lp4 -c debug=1 extensions/lp_call_extension.py -c paramfile=\"data/toy_model_constraint.opt\" extensions/efm_checker.py -c efmcheckfile=\"data/toy_model_biofuel.efmc\" data/toy_model_biofuel.lp4 -n 0 > output_biofuel_efms_lp_constraint.txt
#$clingoLP $clingoLPParams $clingoParams $minSupport extensions/debug_propagator.py -c debug=2 -c trace=1 solveLP.lp4 extensions/lp_call_extension.py -c paramfile=\"data/toy_model_constraint.opt\" extensions/efm_checker.py -c efmcheckfile=\"data/toy_model_biofuel.efmc\" data/toy_model_biofuel.lp4 -n 0 > output_biofuel_efms_lp_constraint.txt
echo "> written in output_biofuel_efms_lp_constraint.txt"
