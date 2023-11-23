#const mcscheckfile=nofile.
#const mcssavefile=nofile.

#script (python)

import itertools
import pickle
import random
import re
import sys
import warnings
import numpy as np
import cobra

from extensions.clingoLP_extension import clingoLPExtension
from extensions.support_propagator import SupportPropagator
import extensions.parsehelper as ph

DEFAULT_DEBUG_VALUE = False # False # debug for printing literals
SHOW_SUPERSETS = True # debug for printing false supersets
WANTED_REACTIONS_SIGNATURE = "wanted" # asp atom signature for wanted reactions

def knock_out(model, lr):
    model = model.copy()
    for r in lr:
        lb = False
        r = r[:] # remove quotes and ending bracket
        if r.startswith('mcs_'):
            r = r[4:]
        if r.startswith('R_'):
            r = r[2:]
        if r.endswith('_rev'):
            r = r[:-4]
            lb = True
        gr = model.reactions.get_by_id(r)
        if lb:
            gr.lower_bound = 0
        else:
            gr.upper_bound = 0
    return model

def sol_after_ko(lr, model):
    try:
        model_c = knock_out(model, lr)
        opt = model_c.optimize()
        return opt.objective_value if opt.status != 'infeasible' else 0.0
    except KeyError as e:
        warnings.warn("Reaction name not found in SBML; MCSChecker likely not invoked with correct SBML file.")
        return 0.0

def nm1_combinations(orig_lr, model, saves):
    all_lr = {}
    for lr in itertools.combinations(orig_lr, len(orig_lr)-1):
        fs = frozenset(lr)
        if fs not in saves:
            saves[fs] = sol_after_ko(lr, model)
        all_lr[fs] = saves[fs]
        if all_lr[fs] == 0.0:
            saves[frozenset(orig_lr)] = 0.0
            return 0.0, all_lr, fs
    orig = sol_after_ko(orig_lr, model)
    saves[frozenset(orig_lr)] = orig
    return orig, all_lr, None
        
def get_reaction_name(word):
    """
    Gets the reaction name from a formatted string in argument
    Returns the regex matching string if search is successful

    Params:
        word: string, clingoLP theory atom flux(v)

    Returns:
        match: a regex match or None if search failed
    """
    res = re.compile(r'flux\(("?.*)"?\)').search(word)
    return None if res is None else res.group(1)    
    
class MCSChecker(clingoLPExtension):

    def __init__(self, mcs_checker_file, stoichreacs, compressed=True, save_file=None):
        self.model = cobra.io.read_sbml_model(mcs_checker_file)
        
        self.total_solutions = self.mcs_tested_solutions = 0
        self.sup_mcs_solutions = self.sub_mcs_solutions = 0
        
        self.stoichreacs = stoichreacs
        assert(self.stoichreacs)

        self.full_last_support = []
        self.last_support = []
        self.literals = {}
        self.saved_sols = dict()
        self.invalid_supports = []
        self.last_state = []
    
        self.save_file = save_file
        self.load_state()

    def filter_support(self, support):
        return [r for r in support if (not r in self.stoichreacs) and (not r.endswith('tgt')) and (not r.endswith('irr'))]
        
    def load_state(self):
        if self.save_file is not None:
            try:
                with open(self.save_file, 'rb') as f:
                    self.saved_sols = pickle.load(f)
                    if not type(self.saved_sols) == dict:
                        warnings.warn('Unable to load save state')
                        self.saved_sols = dict()
            except (OSError, pickle.UnpicklingError):
                warnings.warn('Unable to load save state')
             
    def save_state(self):
        if self.save_file is not None:
            saved_sols = self.saved_sols 
            self.load_state()
            self.saved_sols.update(saved_sols)
            with open(self.save_file, 'wb') as f:
                pickle.dump(self.saved_sols, f)
    
    def get_invalid_support(self, with_inactive=False):
        support = self.invalid_supports.pop()
        active = [self.literals[s] for s in support]
        if with_inactive:
            inactive = [self.literals[r] for r in self.intreacs if r not in support]
            return (active, inactive) 
        return active
    
    def there_are_invalid_supports(self):
        if self.invalid_supports:
            return True
        return False
        
    def is_mcs(self, support):
        support = self.filter_support(support)
        if len(support) < 1:
            warnings.warn('Trying to check if support with no reactions is a MCS')
            return False
        self.last_support = support
        orig, nm1, partial_zero = nm1_combinations(support, self.model, self.saved_sols)
        if (orig == 0.0):
           if partial_zero is not None: 
               self.invalid_supports.append(set(partial_zero))
               return False
           return True
        self.invalid_supports.append(support)
        return False 

    def is_a_mcs(self, lp_state):
        if isinstance(lp_state, str): # only tuple accepted
            return False
        support = [k for k, v in lp_state[1].items() if v > 0]
        support = list(map(get_reaction_name, support))
        support = list(filter(lambda x: x is not None, support))
        support = list(map(ph.smart_remove_quotes, support))
        self.total_solutions += 1
        self.full_last_support = support
        return self.is_mcs(support)
    
    
    def is_valid_partial_solution(self, support):
        support = self.filter_support(support)
        if len(support) < 1:
            warnings.warn('Trying to check if support with no reactions is a valid partial solution')
            return True
        self.last_support = support
        orig, nm1, partial_zero = nm1_combinations(support, self.model, self.saved_sols)
        if (orig == 0.0):
            if partial_zero is not None: 
                self.invalid_supports.append(set(partial_zero))
            else:
                self.invalid_supports.append(support)
            return False
        return True 
    
    def partial_mcs(self, lp_state, wanted):
        if isinstance(lp_state, str): # only tuple accepted
            return False
        support = [k for k, v in lp_state[1].items() if v > 0]
        support = list(map(get_reaction_name, support))
        support = list(filter(lambda x: x is not None, support))
        support = list(map(ph.smart_remove_quotes, support))
        support = list(filter(lambda x: x not in wanted, support))
        self.total_solutions += 1
        self.full_last_support = support
        return self.is_valid_partial_solution(support)
    
    def affected_literals(self, support=[], with_inactive=False):
        if not support:
            support = self.last_support
        active = [self.literals[r] for r in support]
        if with_inactive:
            inactive = [self.literals[r] for r in self.reacs if r not in support]
            return active, inactive
        return active

    def print_stats(self):
        print("Excluded " + str(self.sup_mcs_solutions) + " supersets and " + str(self.sub_mcs_solutions) + " subsets out of " + str(self.total_solutions) + " minimal clingoLP solutions")

    def show_support(self, lits):
        return sorted([self.invlit[abs(lit)] for lit in lits if lit in self.literals.values()])

    @property
    def reacs(self):
        return list(self.literals.keys())
    
    @property
    def intreacs(self):
        return self.filter_support(list(self.literals.keys()))
    
    @property
    def lits(self):
        return list(self.literals.values())

    @property
    def invlit(self):
        return {v: k for k, v in self.literals.items()}


class MCSCheckerExtension(MCSChecker, clingoLPExtension):

    def __init__(self, mcs_checker_file, stoichreacs, wanted, save_file):
        MCSChecker.__init__(self, mcs_checker_file, stoichreacs, wanted, save_file)
        self.debug = DEFAULT_DEBUG_VALUE
        self.wanted = wanted
        print(self.wanted, flush=True)
        
    def init_action(self, init):
        tmp_prop = SupportPropagator()
        tmp_prop.init(init)
        self.literals = tmp_prop.literals()
        assert(self.literals)
        if self.debug:
            tmp_prop.dump_literals()

    def check_consistency_action(self, control, state):
        """ TODO: Implement addition of minimal nogoods wrt constraints instead of non-minimal nogoods """
        if state.current_lit_percentage >= 100:
            if not self.is_a_mcs(state.current_assignment):
                self.sup_mcs_solutions += 1
                assert(self.there_are_invalid_supports()) # if fail here, check filtering method
                active_lits, inactive_lits = self.get_invalid_support(with_inactive=True)
                if SHOW_SUPERSETS:
                    print('Superset filtered', self.show_support(active_lits))
                if not control.add_nogood(active_lits) or not control.propagate():
                    return False
                if not control.add_clause(inactive_lits) or not control.propagate():
                    return False
                warnings.warn("Unable to add nogoods: previously added nogoods were non-minimal wrt constraints")
                return False
            return True
        elif self.wanted and state.current_lit_percentage >= 99:
            if not self.partial_mcs(state.current_assignment, self.wanted):
                self.sub_mcs_solutions += 1
                assert(self.there_are_invalid_supports())
                active_lits, inactive_lits = self.get_invalid_support(with_inactive=True)
                if SHOW_SUPERSETS:
                    print('Partial superset filtered', self.show_support(active_lits))
                if not control.add_nogood(active_lits) or not control.propagate():
                    return False
                if not control.add_clause(inactive_lits) or not control.propagate():
                    return False      
                warnings.warn("Unable to add nogoods: previously added nogoods were non-minimal wrt constraints")
                return False
        return True

    def on_model_action(self, thread_id):
        self.save_state()

    def after_prg_solve_action(self):
        self.print_stats()

    @classmethod
    def from_prg(cls, prg):
        mcs_checker_file_name = str(prg.get_const('mcscheckfile')) 
        if mcs_checker_file_name == 'nofile':
            raise ValueError('MCS Checker invoked without File Name parameter')
        mcs_checker_file = ph.smart_remove_quotes(mcs_checker_file_name)
        cutsets_file_name = str(prg.get_const('mcssavefile')) 
        cutsets_file = ph.smart_remove_quotes(cutsets_file_name) if cutsets_file_name != 'nofile' else None
        stoichreacs = [str(atom.symbol.arguments[0].string) for atom in prg.symbolic_atoms.by_signature("stoichreac", 1)]
        wanted = [str(atom.symbol.arguments[0].string) for atom in prg.symbolic_atoms.by_signature(WANTED_REACTIONS_SIGNATURE, 1)]
        stoichreacs = list(map(ph.smart_remove_quotes, stoichreacs))
        wanted = list(map(ph.smart_remove_quotes, wanted))
        return cls(mcs_checker_file, stoichreacs, wanted, cutsets_file)

#end.
