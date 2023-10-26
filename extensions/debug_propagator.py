#script (python)

from extensions.clingoLP_extension import PropagatorBasedClingoLPExtension
from extensions.support_propagator import SupportPropagator

"""
Most basic inheritance of SupportPropagator from support_propagator.py for debugging purposes.
Toggle debug=2 for full debug of when support literals are set to true and unset.
Let debug=1 if only the correspondance between reaction names and literal numbers is needed.
"""
class DebugPropagator(SupportPropagator):
    def __init__(self, debug_value=1):
        super().__init__()
        self.debug = debug_value

    def init(self, init):
        super().init(init)

"""
Declares the additional extension to be taken into account by clingoLPCplex.py.
"""
class DebugExtension(PropagatorBasedClingoLPExtension):
    def __init__(self, debug_value=1):
        super().__init__(DebugPropagator(debug_value))

    @classmethod
    def from_prg(cls, prg):
        debug_value = int(prg.get_const('debug').number)
        if debug_value < 0 or debug_value > 3:
            raise ValueError('DebugPropagator debug value is incorrect')
        return cls(debug_value)

#end.
