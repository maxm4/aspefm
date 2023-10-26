import sys

from extensions.clingoLP_extension import clingoLPExtension

CST_valid_solution = True
CST_invalid_solution = False

def get(val, default):
    return val if val != None else default

class ListOfExtensions(object):

    def __init__(self, prg, glb):
        self.listOfExtensions = [] # List<clingoLPExtension>
        globals().update(glb) # Adds the extensions imported from clingo script integration
        self.add_extensions(prg) # prg contains the clingo CLI parameters
        self.print_extensions() # print every extensions, enabled by default for output clarity

    def assert_valid_list(self):
        for extension in self.listOfExtensions:
            assert(type(extension) == clingoLPExtension)

    def print_extensions(self):
        for extension in self.listOfExtensions:
            print('Extension', type(extension).__name__, 'loaded')

    """ Defining every action as calling the corresponding method of every extension of the list  """

    def init_action(self, init):
        for ext in self.listOfExtensions:
            ext.init_action(init)

    def decide_action(self, thread_id, assignmt, fallback):
        for ext in self.listOfExtensions:
            fallback = ext.decide_action(thread_id, assignmt, fallback)
        return fallback

    def propagate_before_update_action(self, control, state, changes):
        for ext in self.listOfExtensions:
            if not ext.propagate_before_update_action(control, state, changes):
                return CST_invalid_solution
        return CST_valid_solution

    def check_consistency_action(self, control, state):
        for ext in self.listOfExtensions:
            if not ext.check_consistency_action(control, state):
                return CST_invalid_solution
        return CST_valid_solution

    def propagate_after_check_action(self, control, state):
        for ext in self.listOfExtensions:
            if not ext.propagate_after_check_action(control, state):
                return CST_invalid_solution
        return CST_valid_solution

    def undo_action(self, thread_id, assign, changes):
        for ext in self.listOfExtensions:
            ext.undo_action(thread_id, assign, changes)

    def cplex_init_action(self, cplex):
        for ext in self.listOfExtensions:
            ext.cplex_init_action(cplex)

    def before_prg_solve_action(self, prg):
        for ext in self.listOfExtensions:
            ext.before_prg_solve_action(prg)

    def on_model_action(self, thread_id):
        for ext in self.listOfExtensions:
            ext.on_model_action(thread_id)

    def after_prg_solve_action(self):
        for ext in self.listOfExtensions:
            ext.after_prg_solve_action()

    """ Adding every extension, actions are order sensitive so order is to be considered """

    def add_extensions(self, prg):
        self.add_debug_extension(prg)
        self.add_efm_checker(prg)
        self.add_mcs_checker(prg)
        self.add_cplex_debug(prg)
        self.add_lp_caller(prg)
        self.add_profiler(prg)

    """ Initializing every extension that was found using CLI parameters """

    def add_debug_extension(self, prg):
        try:
            type(DebugExtension)
            debug_ext = DebugExtension.from_prg(prg)
            self.listOfExtensions.append(debug_ext)
        except NameError:
            pass
        
    def add_efm_checker(self, prg):
        try:
            type(EFMCheckerExtension)
            efm_checker = EFMCheckerExtension.from_prg(prg)
            self.listOfExtensions.append(efm_checker)
        except NameError:
            pass

    def add_mcs_checker(self, prg):
        try:
            type(MCSCheckerExtension)
            mcs_checker = MCSCheckerExtension.from_prg(prg)
            self.listOfExtensions.append(mcs_checker)
        except NameError:
            pass

    def add_cplex_debug(self, prg):
        try:
            type(cplexDebugExtension)
            cdb_ext = cplexDebugExtension()
            self.listOfExtensions.append(cdb_ext)
        except NameError:
            pass

    def add_lp_caller(self, prg):
        try:
            type(LPCallExtension)
            lp_caller = LPCallExtension.from_prg(prg)
            self.listOfExtensions.append(lp_caller)
        except NameError:
            pass

    def add_profiler(self, prg):
        try:
            type(ProfilerExtension)
            prf_ext = ProfilerExtension.from_prg(prg)
            self.listOfExtensions.append(prf_ext)
        except NameError:
            pass

