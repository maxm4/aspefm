
import time
import clingo

class SupportPropagator:

    def __init__(self):
        self.cumultime = 0
        self.watch_neglits = False
        self.debug = 0
        self.support_name = "support"
        self.tuple_arity = 1
        self.reactions = {} # shared state, computed by init and unchanged later
        self.states = [] # per thread state, computed and modified by propagator and undo
        self.reaclits = self.reactions.keys()
        self.reacnames = self.reactions.values()
        self.alllits = [] # all watched positive literals
        # reactions is a dict associating support atom literal numbers to reaction names
        # states contains sets of positive watched literals (one per thread)
        
    def dump_literals(self):
        print('Dumping reactions/literals dict...')
        print(self.literals(), flush=True)

    def support(self, state): # get list of active reaction names, from positive watched literals
        return [self.reactions[lit] for lit in state if lit in self.reaclits] 

    def reaction_states(self, state): # get tuples of (reaction, bool, lit) describing active and inactive reactions in state
        return [(self.reactions[abs(lit)], (lit > 0), lit) for lit in state if abs(lit) in self.reaclits] 

    def watched_literals(self, state): # filters list of active watched literals to only those corresponding to support
        return [lit for lit in state if lit in self.reaclits] 

    def literals(self): # gets dict of reactions associated to literals
        return {reaction: lit for lit, reaction in self.reactions.items()}

    def inactive(self, state): # get list of inactive literals, from positive watched literals
        return [-lit for lit in self.alllits if lit not in state]

    def fullstate(self, state): # get complete state for nogood creation
        return [(-1 if lit not in state else 1)*lit for lit in self.alllits] 

    def argname(self, arg):
        if arg.type == clingo.SymbolType.String:
            return arg.string
        if arg.type == clingo.SymbolType.Number:
            return arg.number
        if arg.type == clingo.SymbolType.Function:
            return arg.name
        return arg
    
    def parse_arguments(self, args):
        return "_".join([str(self.argname(x)) for x in args]) 

    def init_action(self, init):
        pass
    
    def propagate_action(self, control, state):
        return True

    def undo_action(self, state):
        pass

    def init(self, init):
        self.states = [ set() for _ in range(init.number_of_threads) ]
        for atom in init.symbolic_atoms:
            if atom.match(self.support_name, self.tuple_arity):
                lit = init.solver_literal(atom.literal)
                init.add_watch(lit)
                if self.watch_neglits: init.add_watch(-lit)
                self.alllits.append(lit)
                self.reactions[lit] = self.parse_arguments(atom.symbol.arguments)
        if self.debug > 0: self.dump_literals()
        self.init_action(init)

    def propagate(self, control, changes):
        start = time.time()
        state = self.states[control.thread_id]
        for lit in changes:
            state.add(lit)
        if self.debug > 1:
            poschanges = [x for x in changes if x > 0]
            if poschanges:
                self.msg_prop_debug(state)
        action = self.propagate_action(control, state)
        end = time.time()
        self.cumultime += end-start
        if self.debug > 2:
            print(action, end=' ')
        return action

    def undo(self, thread_id, assignment, changes):
        start = time.time()
        state = self.states[thread_id]
        for lit in changes:
            state.remove(lit)
        if self.debug > 1:
            poschanges = [x for x in changes if x > 0]
            if poschanges:
                self.msg_undo_debug(poschanges, state)
        self.undo_action(state)
        end = time.time()
        self.cumultime += end-start

    def decide(self, thread_id, assignmt, fallback):
        return fallback

    def msg_prop_debug(self, state):
        print('prop', self.support(state), flush=True)        

    def msg_undo_debug(self, changes, state):
        print('undo', self.support(changes), self.support(state), flush=True)

    def print_assignment(self, thread_id):
        pass  

"""
Main function is in clingoLPCPlex.py
We can only have one main function but we want both propagators to operate at the same time

The idea is the following:

    prop = Propagator(clingo_LP_params)
    add_prop = AdditionalPropagator()
    prg.register_propagator(add_prop)
    prg.register_propagator(prop)
    prg.ground([("base", [])])
    prg.solve(on_model = print_assignment)

Now replaced by the extensions and listOfExtensions system, less time-consuming than additional propagators
"""
