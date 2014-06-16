# -*- coding: utf-8 -*-
"""
HASH input/output core module

Here one can specify and register IO formats for HASH.

More doc to come, basically, each format is a key in this
module's global dict called "IO_REGISTRY", and each value is
a dict containing a key called "module" with the module name,
and optionally "input" and "output" keys specifying the names of
the input and output functions in the module. The default
function names are "input" and "output".


"""
import importlib

### Add I/O formats here! ####################################################
#
IO_REGISTRY = {"OBSPY" : { 'module': "hashpy.io.obspyIO",
                           'input' : "inputOBSPY",
                           'output': "outputOBSPY",
                         },

               "ANTELOPE" : {'module': "hashpy.io.antelopeIO"},
                
               "FPFIT" : {'module': "hashpy.io.fpfitIO"},
              }
##############################################################################


class IOFunction(object):
    """
    Class whose instances are a callable function loaded from a module

    Designed for HASH I/O, but pretty generic. Could break out to a 'util' mod
    """
    _format = None
    _formatmap = dict()
    _registry = IO_REGISTRY
    _module = None
    _fxn = None

    def __init__(self, function, format=None):
        """
        Load a module from a mapping, get a specified function
        """
        if format is not None:
            self._format = format
            if format in self._registry:
                self._formatmap = self._registry[self._format.upper()]
                module_name = self._formatmap.get('module')
            else:
                module_name = self._format 
            self._module = importlib.import_module(module_name) # TODO catch keyerror, importerror
        # Allow mapping to a different function name 
        func_name = self._formatmap.get(function, function)
        if hasattr(self._module, func_name):
            self._fxn = getattr(self._module, func_name)

    def __call__(self, *args, **kwargs):
        return self._fxn(*args, **kwargs)


#
# EXAMPLE used as default if an output format can't be determined
#
def output(hp):
    """
    Simple string line output of best solution
    
    NOTES
    -----
    This is called as a default if hp.output() is called with no format.
    
    Uses the hp._best_quality_index method, from the original HASH code,
    so this is easily modified to a custom quality assessment/output

    """
    x = hp._best_quality_index
    s,d,r = hp.str_avg[x], hp.dip_avg[x], hp.rak_avg[x]
    return 'Solution:{orid} |  STRIKE: {st:0.1f}  DIP: {dp:0.1f}  RAKE: {rk:0.1f}  | Quality:{q}'.format(orid=hp.icusp,
        st=float(s), dp=float(d), rk=float(r), q=hp.qual[x])

