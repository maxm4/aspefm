clingoLP="clingo clingoLP.py"
clingoLPParams=(-c nstrict=0 -c epsilon="(1,2)" -c accuracy=10)
clingoLPParams="${clingoLPParams[@]}"
minSupport="--heuristic Domain --enum-mode domRec"
clingoParams="--stats=2"

echo "== all minimal cut sets of the base network if without target =="
$clingoLP mcsLP.lp4 data/mcs_example_hadicke.lp4 -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_no_target.txt
echo "> written in output_mcs_example_mcs_no_target.txt"

echo "== target: no synthesis of undesired products D and E =="
$clingoLP mcsLP.lp4 data/mcs_example_hadicke.lp4 data/mcs_example_target.txt -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_target.txt
echo "> written in output_mcs_example_mcs_target.txt"

echo "== no synthesis of undesired products D and E; EM1 or EM2 =="
$clingoLP mcsLP.lp4 data/mcs_example_hadicke.lp4 data/mcs_example_target.txt data/mcs_example_constraint_1.txt -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_cstr_1.txt
echo "> written in output_mcs_example_mcs_cstr_1.txt"

echo "== no synthesis of undesired products D and E; EM1 and EM2 =="
$clingoLP mcsLP.lp4 data/mcs_example_hadicke.lp4 data/mcs_example_target.txt data/mcs_example_constraint_2.txt -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_cstr_2.txt
echo "> written in output_mcs_example_mcs_cstr_2.txt"

echo "== same target; Wanted reactions R5 and R10; No MCSChecker =="
$clingoLP mcsLP.lp4 data/mcs_example_hadicke.lp4 data/mcs_example_target.txt data/mcs_example_wanted_constraint.txt -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_wanted_supersets.txt
echo "> written in output_mcs_example_mcs_wanted_supersets.txt"

echo "== same target; Wanted reactions R5 and R10; With MCSChecker =="
$clingoLP mcsLP.lp4 extensions/mcs_checker.py -c mcscheckfile=\"data/mcs_example.xml\" data/mcs_example_hadicke.lp4 data/mcs_example_target.txt data/mcs_example_wanted_constraint.txt -n 0 $clingoLPParams $minSupport $clingoParams > output_mcs_example_mcs_wanted.txt
echo "> written in output_mcs_example_mcs_wanted.txt"




