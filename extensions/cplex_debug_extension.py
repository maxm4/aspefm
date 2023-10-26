
#script (python)
import sys
from extensions.clingoLP_extension import clingoLPExtension

"""
Declares an additional extension for cplex debugging to be taken into account by clingoLPCplex.py.
"""
class cplexDebugExtension(clingoLPExtension):

    def __init__(self, cplexDebug=1, lastlp=None):
        super().__init__()
        self.cplexDebug = cplexDebug
        self.lastlp = lastlp

    def cplex_init_action(self, cplex):
        cplex.__cplex_debug = self.cplexDebug ## VALUE EDIT with influence on clingoLP.py
        self.cplex = cplex.__solver_obj ## WARNING: multithread unsafe
        if self.cplexDebug > 1:
            cplex.__solver_obj.set_log_stream(sys.stderr)
            cplex.__solver_obj.set_error_stream(sys.stderr)
            cplex.__solver_obj.set_warning_stream(sys.stderr)
            cplex.__solver_obj.set_results_stream(sys.stderr)
        return cplex
        
        ## cplex.__solver_obj.parameters.conflict_algorithm = 4 # force conflict refiner to compute irreducibly inconsistent sets
        ## cplex.__solver_obj.parameters.conflict.display = 2 # if logs activated
        ## cplex.__solver_obj.parameters.simplex.tolerances.markowitz = 0.1 # markowitz tolerance

    def after_prg_solve_action(self):
        if not self.cplex or not self.lastlp: 
            return
        self.cplex.write_file(self.lastlp)

#end.
