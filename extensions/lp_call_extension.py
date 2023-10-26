#const paramfile=nofile.

#script (python)

from extensions.clingoLP_extension import PropagatorBasedClingoLPExtension
import extensions.parsehelper as ph
from extensions.support_propagator import SupportPropagator

import numpy as np
import time
import cplex
from cplex.exceptions import CplexSolverError, error_codes

CST_valid_solution = True
CST_invalid_solution = False

defluxer = lambda x: x[len('flux("'):-len('")')]


# TODO: ultimately remove this class
class ThermoClingoLPCplexCaller(): # unused, check class below instead / also using clingoLP functions is not a good idea


    def __init__(self, state):
        self.state = state # from the state object from clingoLP, one can get an instance of cplex
        pass


    def get_conflict(self, state): # cplex conflict refiner if necessary, for addition of clingoLP conflicts
        confl = state.lp.get_confl()
        confl.refine(confl.linear_constraints())
        try:
            confl_cnums = np.array(confl.get())
        except CplexSolverError as err:
            if not err.args[2] == error_codes.CPXERR_NO_CONFLICT:
                raise
            return []
        confl_cnums = np.where(confl_cnums == confl.group_status.member)[0]
        confl_cnums = [self.literals() for i in confl_cnums]
        return confl_cnums
        

    """ Delta G attempt: not tested, not used in rest of the code """
    def compute_dg_lp(self, current_support, state): # attempt to reinitialize clingoLP's cplex caller to use it for our purposes
        state.lp.reset()
        constr_list = []
        self.r_to_dg = {}
        lbs, ubs = [], []
        for i, reaction in enumerate(current_support): # based on support but thermo program needs fluxes
            self.r_to_dg[i] = reaction
            dg_i = {reaction: '1'}
            SjiDG = 0
            for w, k in self.stoch[reaction]:
                dg_i[w] = -k*self.RT # w metab var is log cj/c0
                lbs.append(self.param_file_dict)
                SjiDG += k*self.dgf0[w]
                if w not in lbs and w not in ubs:
                    lbs.append((w, self.mbounds[w][0]))
                    ubs.append((w, self.mbounds[w][1]))
            constr = (dg_i, '=', SjiDG)
            ubs.append((reaction, 0)) # delta rGi <= 0 
            constr_list.append(constr)
        state.lp.add_constr(constr_list)
        state.lp.__solver_obj.variables.set_lower_bounds(lbs)
        state.lp.__solver_obj.variables.set_upper_bounds(ubs)
        state.lp.solve_lp()
        return state.lp.is_sat()


# Main class
class LPCallChecker():

    def __init__(self, param_file):
        """
        LPCallChecker class constructor

        Params:
            param_file (str): input file name, file contains a dict of costs for each reaction
        """
        super().__init__()
        with open(param_file, 'r') as f:
            self.data = eval(f.read())

        self.constraints = self.data
        self.added_nogoods = 0 # incremented in extension
        self.checked_solutions = 0 # incremented in extension
        self.cumultime = 0 # to be incremented in respects_constraints        
    
    " Basic idea where a sum is performed based on the flux value of the partial solution, and tested against a bound"
    " Ultimately, this can be replaced by a full LP call, based on the flux values of the partial solution* "
    " *If there is a real optimization problem to be done, with undetermined variables (ex: thermo concentrations) "
    def respects_constraints(self, fluxes): # can be used for testing if above max enzyme requirements for example
        cumultime = time.time()
        for constraint in self.constraints:
            linsum = sum(coeff * fluxes[reac] for coeff, reac in constraint)
            if not linsum <= 0: # how constraints are defined in the current data
                return False
        return True # only if all constraints respected

    def respects_lp_call(self, state):
        fluxes = self.get_fluxes(state.current_assignment)
        #print(fluxes)
        deb = time.time()
        bool_respected = self.respects_constraints(fluxes)
        end = time.time()
        self.cumultime += (end - deb)
        return bool_respected

    def print_stats(self):
        print("Added nogoods for " + str(self.added_nogoods) + " out of " + str(self.checked_solutions) + " partial solutions that were inadequate")
        print("Total time used by LPCallExtension: " + str(round(self.cumultime, 5)))

    def get_support(self, current_assignment):
        return [defluxer(k) for k, v in current_assignment[1].items() if v != 0]

    def get_fluxes(self, current_assignment):
        return {defluxer(k): v for k, v in current_assignment[1].items()}

    def get_literals(self, support): # return positive literals for support active reactions
        return [self.literals[s] for s in support] # self.literals filled in extension

    def get_clause(self, support):
        active = [self.literals[r] for r in support]
        inactive = [self.literals[r] for r in self.reacs if r not in support]
        return active, inactive

    @property
    def reacs(self):
        return list(self.literals.keys())

# Extension for the list of extensions
class LPCallExtension(LPCallChecker, clingoLPExtension):

    def __init__(self, param_file, debug_value):
        LPCallChecker.__init__(self, param_file)
        self.debug = (debug_value > 0) 

    # TODO: transform SupportPropagator into SupportExtension, integrate prop.init in ClingoLP init 
    def init_action(self, init): 
        tmp_prop = SupportPropagator()
        tmp_prop.init(init)
        self.literals = tmp_prop.literals()
        assert(self.literals)
        tmp_prop.dump_literals()

    def check_consistency_action(self, control, state):
        if not self.respects_lp_call(state):
            self.checked_solutions += 1
            support = self.get_support(state.current_assignment)
            if self.debug: print('lp support', support)
            literals = self.get_literals(support)
            if self.debug: print('lits for clause', literals)
            active, inactive = self.get_clause(support)
            if self.debug: print(inactive)
            # 2 possible options: add nogood on active, or add clause on inactive, when possible
            if not control.add_nogood(active, lock=True) or not control.propagate():
                self.added_nogoods += 1
                return False
            # I believe this two logical propositions should be equivalent, in doubt always use the nogood on active version
            if not control.add_clause(inactive, lock=True) or not control.propagate():
                self.added_nogoods += 1
                return False
            warnings.warn("Potentially abnormal behaviour: unable to add nogoods")
            # To solve this I believe the only solution is to get the core conflict, like is done with the stoichiometry constraint in clingoLP
            # In the case of data/toy_model_constraint.opt the conflicting constraint is simply R2, since R2 - 3 R10 <= 0 and no EFM has this two together, all EFMs containing R2 are eliminated
            # So adding a second check in clingoLP, yes, but we need a way to get core conflicting constraints instead of just conflicting constraints
            # Using the stoichiometry constraints in clingoLP code and the cplex refiner would be a good way of achieving this
            # The problem being that this does not stop propagation, since the solver does not know that R2 is the sole responsible
            # Propagation continues and cplex keeps calculating the same LP solution, until the constraint (R2 = 0) is finally propagated
            return False
        return True

    def after_prg_solve_action(self):
        self.print_stats()

    @classmethod
    def from_prg(cls, prg):
        param_file_name = str(prg.get_const('paramfile')) 
        debug_value = int(prg.get_const('debug').number)
        if param_file_name == 'nofile':
            raise ValueError('LP Call extension invoked without File Name parameter')
        param_file_name = ph.smart_remove_quotes(param_file_name)
        return cls(param_file_name, debug_value)

#end.
