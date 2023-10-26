#const efmcheckfile=nofile.

#script (python)

import pickle
import re
import numpy as np

from extensions.clingoLP_extension import clingoLPExtension
import extensions.parsehelper as ph
from extensions.support_propagator import SupportPropagator

HEUR_REACTIONS = 70 # percentage threshold
DEFAULT_DEBUG_VALUE = False # debug for printing literals
SHOW_SUPERSETS = False # False # True # debug for printing false supersets

def is_kernel_dimension_one(matrix):
    """
    Checks if the kernel of a matrix is of dimension one
    Uses SVD to calculate the rank of the matrix

    Params:
        matrix: stoichiometric matrix

    Returns:
        boolean: True if the kernel is of dimension one, else False
    """
    # utilize SVD instead of Gaussian Elimination (done in Klamt, 2005)
    rank = np.linalg.matrix_rank(matrix)
    # print("forme matrice", matrix.shape, "rang", rank)
    if (rank <= matrix.shape[1] - 1):  # nb columns aka nb reactions
        if (rank == matrix.shape[1] - 1):
            return True # efm
        return None # subset of efm
    return False # superset of efm


def possible_efm_according_to_kernel(matrix, unbound):
    """
    Checks if an efm is possible according to the kernel rank

    Params:
        matrix: stoichiometric matrix

    Returns:
        boolean: True if the kernel is of dimension one, else False
    """
    # utilize SVD instead of Gaussian Elimination (done in Klamt, 2005)
    rank = np.linalg.matrix_rank(matrix)
    # print("forme matrice", matrix.shape, "rang", rank)
    if (rank <= matrix.shape[1] - 1):  # nb columns aka nb reactions
        if (rank == matrix.shape[1] - 1):
            return True # efm
        if (unbound > 0) and (rank >= (matrix.shape[1] - 1) - unbound):
            print("Not actually filtered", rank, matrix.shape[1] - 1, rank - (matrix.shape[1] - 1), unbound)
            return True # subset of efm
    if SHOW_SUPERSETS:
        print("MCFM of level", rank, matrix.shape[1] - 1, rank - (matrix.shape[1] - 1), unbound)
    return False # superset of efm


def support_as_boolean(support, reacidx):
    """
    Converts support to a boolean table matching the network file matrix

    Params:
        support: set of string names, support of the reaction, from clingoLP output
        reacidx: hash map associating reaction names to indices, from network file

    Returns:
        booltable: table of booleans:
                   True if the reaction is in the support, else False
    """
    booltable = []
    itersup = support.copy()
    for val in reacidx: # for each reaction name
        match = val in itersup # boolean
        booltable.append(match)
        if match:
            itersup.remove(val)
    return booltable


def is_efm(support, reacidx, matrix):
    """
    Checks if a support reaction names set is an elementary flux mode

    Params:
        fluxv: a flux vector returned by clingo[LP], as pandas dataframe
        reacidx: Python Dict of reaction names corresponding to matrix indices
        matrix: original matrix, to check flux vectors minimality

    Returns:
        status: boolean, true if it is an efm else false
    """
    support = handle_revs(support)
    boolsupp = support_as_boolean(support, reacidx)
    submatrix = matrix[:,boolsupp]
    return is_kernel_dimension_one(submatrix)

def can_be_efm(support, reacidx, matrix, unbound):
    """
    Checks if a support reaction names set can be an elementary flux mode

    Params:
        fluxv: a flux vector returned by clingo[LP], as pandas dataframe
        reacidx: Python Dict of reaction names corresponding to matrix indices
        matrix: original matrix, to check flux vectors minimality
        unbound: number of unbound literals, ie. still unknown matrix dimensions

    Returns:
        status: boolean, true if it is an efm else false
    """
    support = handle_revs(support)
    boolsupp = support_as_boolean(support, reacidx)
    submatrix = matrix[:,boolsupp]
    return possible_efm_according_to_kernel(submatrix, unbound)


def handle_revs(support):
    """
    Removes 'rev' suffixes from the support names

    Params:
        support: list of string names, support of the reaction, from clingoLP 

    Returns:
        support: list of string names, support of the reaction, from clingoLP
    """
    nsupport = []
    for reac in support:
        if ph.is_reversible(reac):
            reac_fwd = ph.remove_rev(reac)
            nsupport.append(reac_fwd)
        else:
            nsupport.append(reac)
    return nsupport


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


class EFMChecker():

    def __init__(self, efm_checker_file):
        
        with open(efm_checker_file, 'rb') as f:
            config = pickle.load(f)

        self.matrix = config.matrix
        self.rindex = config.rindex
        self.neighbours = config.neighbours
        self.neighbours = {k: list(v) for k, v in self.neighbours.items()}
        
        self.nb_reactions = self.matrix.shape[1]
        self.heur_reactions = self.nb_reactions * HEUR_REACTIONS / 100
        self.total_solutions = self.efm_tested_solutions = 0
        self.sup_efm_solutions = 0

        self.last_support = []
        self.checked_supports = set()
        self.literals = {}

    def next_neighbour(self, last_lit, fallback, assignmt):
        sign = 1 if last_lit > 0 else -1
        last_lit = abs(last_lit)
        if last_lit not in self.invlit:
            return fallback
        last_rname = ph.smart_remove_revs(self.invlit[last_lit])
        if self.neighbours[last_rname]:
            for neighbour in self.neighbours[last_rname]:
                if neighbour in self.literals and assignmt.value(self.literals[neighbour]) is None:
                    return sign*self.literals[neighbour]
        return fallback

    def check_efm(self, support):
        if str(support) in self.checked_supports:
            return None
        self.total_solutions += 1
        if len(support) < self.heur_reactions and len(support) > 1:
            self.checked_supports.add(str(support))
            support = list(map(ph.smart_remove_quotes, support))
            self.last_support = support
            support = set(list(map(ph.smart_remove_revs, support)))
            self.efm_tested_solutions += 1
            return is_efm(support, self.rindex, self.matrix)
        return None

    def kamp_to_bal(self, support): 
        for x in support: # converts solution from kamp formalism to bal formalism
            if x.endswith('_rev'):
                irr_isozyme = ph.remove_rev(x) + '_irr'
                if irr_isozyme in self.rindex:
                    support.append(irr_isozyme)
        return support

    def is_an_efm(self, lp_state):
        if isinstance(lp_state, str): # only tuple accepted
            return False
        support = [k for k, v in lp_state[1].items() if v > 0]
        support = list(map(get_reaction_name, support))
        support = list(map(ph.smart_remove_quotes, support))
        support = self.kamp_to_bal(support)
        self.total_solutions += 1
        self.last_support = support
        return is_efm(support, self.rindex, self.matrix)
    
    def can_be_an_efm(self, lp_state, unbound_lits):
        if isinstance(lp_state, str): # only tuple accepted
            return False
        support = [k for k, v in lp_state[1].items() if v > 0]
        support = list(map(get_reaction_name, support))
        support = list(map(ph.smart_remove_quotes, support))
        support = self.kamp_to_bal(support)
        self.total_solutions += 1
        self.last_support = support
        return can_be_efm(support, self.rindex, self.matrix, unbound_lits)

    def is_new_efm(self, lp_state):
        if isinstance(lp_state, str): # only tuple accepted
            return False
        support = [k for k, v in lp_state[1].items() if v > 0] # cplex returns tuple
        if str(support) in self.checked_supports: # already checked
            return False
        self.total_solutions += 1
        if len(support) < self.heur_reactions and len(support) > 1:
            self.checked_supports.add(str(support))
            support = list(map(get_reaction_name, support))
            support = list(map(ph.smart_remove_quotes, support))
            self.last_support = support
            support = self.kamp_to_bal(support)
            support = set(list(map(ph.smart_remove_revs, support)))
            self.efm_tested_solutions += 1
            return is_efm(support, self.rindex, self.matrix)
        return False

    def affected_literals(self, support=[]):
        if not support:
            support = self.last_support
        active = [self.literals[r] for r in self.reacs if r in support]
        inactive = [-self.literals[r] for r in self.reacs if r not in support]
        return active, inactive

    def print_stats(self):
        print("Excluded " + str(self.sup_efm_solutions) + " solutions out of " + str(self.total_solutions) + " minimal clingoLP solutions")

    def show_support(self, lits):
        return sorted([self.invlit[abs(lit)] for lit in lits if lit in self.literals.values()])

    @property
    def reacs(self):
        return list(self.literals.keys())

    @property
    def invlit(self):
        return {v: k for k, v in self.literals.items()}


class EFMCheckerExtension(EFMChecker, clingoLPExtension):

    def __init__(self, efm_checker_file):
        EFMChecker.__init__(self, efm_checker_file)
        self.debug = DEFAULT_DEBUG_VALUE

    def init_action(self, init):
        tmp_prop = SupportPropagator()
        tmp_prop.init(init)
        self.literals = tmp_prop.literals()
        assert(self.literals)
        if self.debug:
            tmp_prop.dump_literals()

    def check_consistency_action(self, control, state):
        if state.current_lit_percentage >= 100:
            unbound_lits = state.total_lits - state.lits_current
            if not self.can_be_an_efm(state.current_assignment, unbound_lits):
                self.sup_efm_solutions += 1
                active_lits, inactive_lits = self.affected_literals()
                if SHOW_SUPERSETS:
                    print('Superset of efm filtered', self.show_support(active_lits))
                if not unbound_lits: 
                    if not control.add_nogood(active_lits) or not control.add_clause(inactive_lits) or not control.propagate():
                        return False   
                    assert(False)
                else:
                    return False
        return True  

    def after_prg_solve_action(self):
        self.print_stats()

    @classmethod
    def from_prg(cls, prg):
        efm_checker_file_name = str(prg.get_const('efmcheckfile')) 
        if efm_checker_file_name == 'nofile':
            raise ValueError('EFM Checker invoked without File Name parameter')
        efm_checker_file = ph.smart_remove_quotes(efm_checker_file_name)
        return cls(efm_checker_file)

#end.
