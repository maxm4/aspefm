#const profilerfile=nofile.
#script (python)

from extensions.clingoLP_extension import clingoLPExtension
import cProfile, pstats, io
"""
Declares an additional extension for the profiler to be taken into account by clingoLPCplex.py.
"""
class ProfilerExtension(clingoLPExtension):

    def __init__(self, dumpfile=None):
        super().__init__()
        self.dumpfile = dumpfile

    def before_prg_solve_action(self, prg):
        self.profile = cProfile.Profile()
        self.profile.enable()

    def after_prg_solve_action(self):
        self.profile.disable()
        if self.dumpfile:
            self.profile.dump_stats(self.dumpfile)
        else:
            self.profile.print_stats()

    @classmethod
    def from_prg(cls, prg):
        param_file_name = str(prg.get_const('profilerfile'))
        param_file_name = ph.smart_remove_quotes(param_file_name)
        if param_file_name == 'nofile':
            param_file_name = None
        return cls(param_file_name)

#end.
