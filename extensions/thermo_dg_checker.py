#const sbml=nofile.
#const concentrations=nofile.
#const deltaG=nofile.

#script (python)

#clingo [...] thermo_dg_checker.py -c debug=1 -c sbml=data.sbml -c concentrations=data.cc -c deltaG=data.dG

from types import SimpleNamespace
from extensions.clingoLP_extension import clingoLPExtension
import extensions.parsehelper as ph
from extensions.support_propagator import SupportPropagator

import numpy as np
import time, warnings

defluxer = lambda x: x[len('flux("'):-len('")')]

class ThermoDGChecker():
    
    def __init__(self, sbml, concentrations, deltaG): # MODIFY
        """
        ThermoDGChecker class constructor

        Params:
            sbml (str): input file name
            concentrations (str): input file name
            deltaG (str): input file name
        """
        with open(sbml, 'r') as f: # SBML
            self.sbml = eval(f.read())  # cobra.io.read_sbml_model
            
        # create model, get metabolites, reactions, stoichiometry
        # reactions sbml DOIVENT AVOIR LE MEME NOM que reactions clingo
        
        # a remplir
        self.metabolites = []
        self.reactions = []
        self.stoichiometry = []
            
        with open(concentrations, 'r') as f: # Concentrations
            self.concentrations = eval(f.read()) 
            
        with open(deltaG, 'r') as f: # DeltaG
            self.deltaG = eval(f.read())  # obtained from Equilibrator
            
        # lecture sbml
        # lecture concentrations min/max
        # lecture deltag formation from equilibrator
        
        self.literals = {} # dict of literals, filled in init_action  
        
        # statistics of checker # DON'T TOUCH
        self.added_nogoods = 0 # incremented in extension
        self.checked_solutions = 0 # incremented in extension
        self.cumultime = 0 # to be incremented in respects_constraints
        
        
    def init_action(self, init): # DON'T TOUCH
        """
        Init action from clingoLPExtension/listOfExtensions.init_action

        Params:
            init: PropagatorInit from clingo/clingoLP
        """
        # Propagator is needed for getting literals
        tmp_prop = SupportPropagator()
        tmp_prop.init(init)
        self.literals = tmp_prop.literals()
        assert(self.literals)
        tmp_prop.dump_literals()   
        
        """
        appelé quand clingoLP.init
        self.literals depend de chaque execution de clingolp, on peut recuperer au moment clingoLP.init
        support EFM = reactions actives
        thermo Checker me dit : support EFM pas bien
        donc on ne veut plus jamais support EFM
        
        ex: R1, R3, on veut contrainte [R1, R3] (ensemble) non autorisé
        ---> nogood
        
        {15: 'R1', 12: 'R3'}
        {'R1': 15, 'R3': 12}
        R1 avec R3 pas autorisé = nogood
        [literal[R1],  literal[R3]]
        15, 12
        # add nogood prend liste en argument
        self.add_nogood([15, 12])
        self.add_nogood([self.literals[r] for r in support])
        """
        
    def construct_constraints(self, efm): # MODIFY
        """
        Construct constraints

        Params:
            efm (dict): reaction name: 1 if active, else 0 (get_support_dict)
        """
        # efm : reactions DOIVENT etre les memes que self.reactions récupérés du SBML        
        print(efm)
        #print(self.metabolites)
        for r in self.reactions:
            if efm[r]  == 1:
                pass
        # print
        return {} # model with constraints


    def call_model_cplex(self, efm): # MODIFY
        """
        Calls CPLEX model to determine if efm is thermodynamically feasible or not

        Params:
            efm (dict): reaction name: 1 if active, else 0 (get_support_dict)
            
        Returns: 
            True if feasible else False
        """       
        self.model = self.construct_constraints(efm)
        sol = SimpleNamespace(status='OPTIMAL')
        #sol = self.model.solve() # type: ignore # docplex.Model
        print(sol)

        if sol.status == 'OPTIMAL':
            return True
        else:
            return False
        

    def respects_thermodynamics(self, state):
        """
        Gets EFM from propagator state
        Calls CPLEX model to determine if efm is thermodynamically feasible or not
        
        Params:
            state (ClingoLP.State): propagator state, to get fluxes and active reactions
            
        Returns: 
            True if feasible else False
        """     
        fluxes = self.get_support_dict(state.current_assignment)
        #print(fluxes)
        deb = time.time()
        bool_respected = self.call_model_cplex(fluxes)
        end = time.time()
        self.cumultime += (end - deb)
        return bool_respected

    """ Utilitary functions to get fluxes or support of EFMs """
    
    def get_support_dict(self, current_assignment):
        return {defluxer(k): 1 if v != 0 else 0 for k, v in current_assignment[1].items()}

    def get_support_list(self, current_assignment):
        return [defluxer(k) for k, v in current_assignment[1].items() if v != 0]

    def get_fluxes(self, current_assignment):
        return {defluxer(k): v for k, v in current_assignment[1].items()}
    
    """ Utilitary functions to get literals from reaction names using propagator literals """

    def get_literals(self, support): # return positive literals for support active reactions
        return [self.literals[s] for s in support] # self.literals filled in extension

    def get_active_list(self, support): # active reactions, not null in efm
        active = [self.literals[r] for r in support]
        return active # list on which we want a nogood
    
    def get_inactive_list(self, support): # inactive reactions, reactions not in support
        inactive = [self.literals[r] for r in self.reacs if r not in support]
        return inactive # list on which we want a clause
    
    @property
    def reacs(self):
        return list(self.literals.keys())
    
    """ 
        Clause: R1 or R2 or R3
        Littéraux: R1, R2, R3 (Ri, ... Rn)
        Opérateur: or
        
        Nogood: (not R1) or (not R2) or (not R3)
        Equivalent: not (R1 and R2 and R3)
        (via De Morgan); on ne veut ni R1, ni R2, ni R3
        Opérator: or
    """
    
    """ Print statistics of Checker"""
    def print_stats(self):
        print("Added nogoods for " + str(self.added_nogoods) + " out of " + str(self.checked_solutions) + " partial solutions checked")
        print("Total time used by ThermoDGChecker: " + str(round(self.cumultime, 5)) + " s")

# Extension for the list of extensions
class ThermoDGCheckerExtension(ThermoDGChecker, clingoLPExtension):

    def __init__(self, sbml, concentrations, deltaG, debug_value):
        ThermoDGChecker.__init__(self, sbml, concentrations, deltaG)
        self.debug = (debug_value > 0) 

    def init_action(self, init): 
        ThermoDGChecker.init_action(self, init)

    def check_consistency_action(self, control, state):
        """
        Check consistency action from clingoLPExtension/listOfExtensions.check_consistency_action
        Is called when consistency is checked after a literal is propagated

        Params:
            control: PropagatorControl is a structure defined in Clingo API that can add nogoods
            state: ClingoLP.State is a structure defined in ClingoLP to store results of propagation and of Cplex
        """
        #if state.current_lit_percentage >= 100: # uncomment this to check total solutions instead of partial
        self.checked_solutions += 1
        if not self.respects_thermodynamics(state):
            support = self.get_support_list(state.current_assignment)
            if self.debug: 
                print('lp support', support)
            active_literals = self.get_active_list(support)
            if self.debug: 
                print('lits for clause', active_literals)
            inactive_literals = self.get_inactive_list(support)
            # DO NOT TOUCH
            if not control.add_nogood(active_literals, lock=True) or not control.propagate():
                self.added_nogoods += 1
                return False
            elif not control.add_clause(inactive_literals, lock=True) or not control.propagate():
                warnings.warn("Adding clause of inactive literals instead of nogoods")
                self.added_nogoods += 1
                return False
            warnings.warn("Potentially abnormal behaviour: unable to add nogoods")
            return False
        return True

    def after_prg_solve_action(self):
        """
        After prg solve action from clingoLPExtension/listOfExtensions.after_prg_solve_action
        Is called after ClingoLP.prg solve is finished
        """
        self.print_stats()
        

    @classmethod
    def from_prg(cls, prg):
        """
        Class method: is called like this: ThermoDGChecker.from_prg(prg)
        
        Params:
            prg: Program arguments, including constants (e.g. paramfile, debug, sbml, ...)
        
        Returns: 
            cls: Class instance using __init__ constructor
        """        
        debug_value = int(prg.get_const('debug').number)
        sbml = ph.smart_remove_quotes(str(prg.get_const('sbml'))) 
        if sbml == 'nofile':
            raise ValueError('ThermoDG Checker invoked without SBML parameter')
        concentrations = ph.smart_remove_quotes(str(prg.get_const('concentrations')))
        if concentrations == 'nofile':
            raise ValueError('ThermoDG Checker invoked without Concentrations parameter')
        deltaG = ph.smart_remove_quotes(str(prg.get_const('deltaG')))
        if deltaG == 'nofile':
            raise ValueError('ThermoDG Checker invoked without DeltaG parameter')
        return cls(sbml, concentrations, deltaG, debug_value)


#end.
