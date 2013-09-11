# -*- coding: utf-8 -*-
#
#  dbhash.py
#
# functions which run HASH using Antelope and hashpy
#
#from obspy.core import UTCDateTime
from hashpy.hashpype import HashPype, HashError
from hashpy.io.antelopeIO import ( load_pf, readANTELOPE, eventfocalmech2db, get_first_motions, RowPointerDict )

import warnings
warnings.filterwarnings("ignore")

def dbhash_run(args):
    """
    Perform a HASH run using Database input and command line args
    """
    hp = HashPype()
    # Load settings data from a pf file...
    if args.pf:
        load_pf(hp,pffile=args.pf)
    else:
        load_pf(hp)
    # Grab data from the db...
    hp.input(args.dbin, format="ANTELOPE", orid=args.orid)
    # Run and catch errors from the minimum requirements checks
    try:
        hp.driver2(check_for_maximum_gap_size=False)
    except HashError as e:
        print "Failed! " + e.message
    except:
        raise
    # For interactive scripting and debugging:
    return hp
    


def main():
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument("dbin",   help="Input database")
    parser.add_argument("dbout",  help="Output database", nargs='?')
    parser.add_argument("-p", "--plot", help="Plot result", action='store_true')
    parser.add_argument("-l", "--loc", help="dbloc2 mode", action='store_true')
    parser.add_argument("--pf",   help="Parameter file")
    group = parser.add_mutually_exclusive_group() #required=True)
    group.add_argument("--evid", help="Event ID", type=int)
    group.add_argument("--orid", help="Origin ID", type=int)
    args = parser.parse_args()
    
    if args.loc:
        hash_prog = dbhash_loc2
    else:
        hash_prog = dbhash_run
    
    hp = hash_prog(args)
    
    # Grab waveform data and launch plotter or spit out solution...
    if args.plot:
        from hashpy.plotting.focalmechplotter import FocalMechPlotter
        # little script to get the FM pick waveform data for plotting
        #adb = get_first_motions(args.dbin, orid=args.orid)
        #t0 = UTCDateTime(RowPointerDict(adb, record=0)['arrival.time']) - 10
        #t1 = UTCDateTime(RowPointerDict(adb, record=adb.nrecs()-1 )['arrival.time']) + 10
        #st = readANTELOPE(adb, starttime=t0, endtime=t1)
        #adb.close()
        #if args.dbout:
        #    savedb = args.dbout
        #else:
        #    savedb = args.dbin       
        ev = hp.output(format="OBSPY")
        p = FocalMechPlotter(ev)
    else:
        # quick orid/strike/dip/rake line
        print hp.output()

    if args.dbout:
        db = hp.output(format="ANTELOPE", dbout=args.dbout)
    
    return hp


if __name__ == '__main__':
    hp = main()    


