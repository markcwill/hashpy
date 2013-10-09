#
"""
HASH input/output core module

Here one can specify and register IO formats for HASH.

More doc to come, basically, each format is a key in this
module's global dict called "IO_REGISTRY", and each value is
a dict containing a key called "module" with the module name,
and optionally "in" and "out" keys specifying the names of
the input and output functions in the module. The default
function names are "input" and "output".


"""
import importlib

### Add I/O formats here! ####################################################
#
IO_REGISTRY = { "OBSPY" : { "module" : "obspyIO" ,
                            "in"     : "inputOBSPY",
                            "out"    : "outputOBSPY",
                          },

                "ANTELOPE" : { "module" : "antelopeIO" },
                
                "FPFIT" : { "module" : 'fpfitIO' },
              }
##############################################################################

#
# EXAMPLE used as default if an output format can't be determined
#
def outputSTRING(hp):
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


class Inputter(object):
    """
    Class whose instances are input/load/read functions for HASH

    """
    
    __input_fxn = None

    @property
    def _input(self):
        return self.__input_fxn

    @_input.setter
    def _input(self, input_function):
        self.__input_fxn = input_function


    def __init__(self, format=None):
        """
        Get the input function and return an inputter that calls it
        """
        if format is not None and format in IO_REGISTRY:
            io_format = IO_REGISTRY[format.upper()]
            io_module = importlib.import_module("hashpy.io." + io_format["module"]) # TODO catch keyerror, importerror
            input_fxn_name = io_format.get("in", "input")
            self._input = getattr(io_module, input_fxn_name)
        else:
            raise NotImplementedError("Can't determine format, must explicity state")

    def __call__(self, *args, **kwargs):
        return self._input(*args, **kwargs)


class Outputter(object):
    """
    Class whose instances are output/write functions for HASH

    """
    
    __output_fxn = None

    @property
    def _output(self):
        return self.__output_fxn

    @_output.setter
    def _output(self, output_function):
        self.__output_fxn = output_function


    def __init__(self, format=None):
        """
        Get the input function and return an inputter that calls it
        
        """
        if format is not None and format in IO_REGISTRY:
            io_format = IO_REGISTRY[format.upper()]
            io_module = importlib.import_module("hashpy.io." + io_format["module"]) # TODO catch keyerror, importerror
            output_fxn_name = io_format.get("out", "output")
            self._output = getattr(io_module, output_fxn_name)
        else:
            self._output = outputSTRING
    
    def __call__(self, *args, **kwargs):
        return self._output(*args, **kwargs)



