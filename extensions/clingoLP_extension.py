

from extensions.support_propagator import SupportPropagator

CST_valid_solution = True
CST_invalid_solution = False

class clingoLPExtension(object):

    def __init__(self):
        pass

    def init_action(self, init):
        pass

    def decide_action(self, thread_id, assignmt, fallback):
        return fallback

    def propagate_before_update_action(self, control, state, changes):
        return CST_valid_solution

    def check_consistency_action(self, control, state):
        return CST_valid_solution

    def propagate_after_check_action(self, control, state):
        return CST_valid_solution

    def undo_action(self, thread_id, assign, changes):
        pass

    def cplex_init_action(self, cplex):
        pass

    def before_prg_solve_action(self, prg):
        pass

    def on_model_action(self, thread_id):
        pass

    def after_prg_solve_action(self):
        pass


class PropagatorBasedClingoLPExtension(clingoLPExtension):

    def __init__(self, prop: SupportPropagator):
        self.prop = prop

    def init_action(self, init):
        self.prop.init(init)

    def decide_action(self, thread_id, assignmt, fallback):
        return self.prop.decide(thread_id, assignmt, fallback)

    def propagate_before_update_action(self, control, state, changes):
        if not self.prop.propagate(control, changes):
            return CST_invalid_solution
        return CST_valid_solution

    def undo_action(self, thread_id, assign, changes):
        return self.prop.undo(thread_id, assign, changes)
    
    def on_model_action(self, thread_id):
        self.prop.print_assignment(thread_id)
