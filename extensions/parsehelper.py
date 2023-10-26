# -*- coding: utf-8 -*-

"""
parsehelper.py

Helper module for mparser
Contains utility functions related to parsing

"""

remove_firstlast = lambda w : w[1:-1] # returns string with no first and last char
in_quotes = lambda w : w.startswith('"') and w.endswith('"')
add_quotes = lambda w : '"{}"'.format(w) # returns string with added quotes
smart_remove_quotes = lambda w : remove_firstlast(w) if in_quotes(w) else w
smart_add_quotes = lambda w : add_quotes(w) if not in_quotes(w) else w
is_reversible = lambda w: smart_remove_quotes(w).endswith('_rev')
smart_remove_revs = lambda w : w[:-4] if w.endswith('_rev') else w
remove_rev = lambda w : w[:-4]
add_rev = lambda w: w + '_rev'

in_rprefix = lambda w : w.startswith('R_')
in_gprefix = lambda w : w.startswith('G_')
in_mprefix = lambda w : w.startswith('M_')
add_rprefix = lambda w : 'R_{}'.format(w)
add_gprefix = lambda w : 'G_{}'.format(w)
add_mprefix = lambda w : 'M_{}'.format(w)
smart_add_rprefix = lambda w : add_rprefix(w) if not in_rprefix(w) else w
smart_add_gprefix = lambda w : add_gprefix(w) if not in_gprefix(w) else w
smart_add_mprefix = lambda w : add_mprefix(w) if not in_mprefix(w) else w
remove_prefix = lambda w : w[2:]
smart_remove_prefix = lambda w : w[2:] if (in_rprefix(w) or in_gprefix(w) or in_mprefix(w)) else w
remove_cpartm = lambda w : w[:-2]
smart_remove_cpartm = lambda w : w[:-2] if (w[-2] == '_' ) else w

def iswfloat(element):
    """
    Returns if the given string can be converted to a float
    
    Params:
          element: float text
    Retruns:
          isfloat: boolean
    """
    try:
        float(element)
        return True
    except ValueError:
        return False


def iswint(element):
    """
    Returns if the given string can be converted to an int
    
    Params:
          element: float text
    Retruns:
          isfloat: boolean
    """
    try:
        int(element)
        return True
    except ValueError:
        return False


def indexl(lst, word):
    """
    Safe index for Python lists
    
    Params:
        lst: list, array
        word: researched word
        
    Returns:
        int: the position of the substring, or -1 if not found
    """    
    try:
        return lst.index(word)
    except ValueError:
        return -1   
