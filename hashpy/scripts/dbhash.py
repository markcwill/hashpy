# -*- coding: utf-8 -*-
"""
#  dbhash.py
#
# functions which run HASH using Antelope and hashpy
#
# optionally, can plot by passing ObsPy output to a plotter
"""
#from obspy.core import UTCDateTime
from hashpy.hashpype import HashPype, HashError
from hashpy.io.antelopeIO import ( load_pf, readANTELOPE, eventfocalmech2db, get_first_motions, RowPointerDict )

def dbhash_run(args):
    """
    Perform a HASH run using Database input and command line args
    
    Input
    -----
    args : Namespace of command line args
    
    Returns : hp : a hashpy.HashPype object containing solutions

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
    
    return hp
    
def main():
    from argparse import ArgumentParser
    
    # Get command line args
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
    
    if args.dbout:
        savedb = args.dbout
    else:
        savedb = args.dbin       
    

    # Now that we have a save location from command line args,
    # make a function to save to that database. The plotter is I/O
    # agnostic, it will accept a function to save anything anyhow anywhichway
    #
    def save_to_db(fmplotter, dbname=savedb, *args, **kwargs):
        focal_mech = fmplotter.event.focal_mechanisms[fmplotter._fm_index]
        if focal_mech is not fmplotter.event.preferred_focal_mechanism():
            fmplotter.event.preferred_focal_mechanism_id = focal_mech.resource_id.resource_id
        
        eventfocalmech2db(event=fmplotter.event, database=dbname)
    

    # Run HASH (special parsing mode for dbloc2)
    if args.loc:
        pass
        # alter args b/c dbloc2 passes a db and a row number
    
    hp = dbhash_run(args)
    

    # Launch plotter or spit out solution
    if args.plot:
        from hashpy.plotting.focalmechplotter import FocalMechPlotter
        ev = hp.output(format="OBSPY")
        p = FocalMechPlotter(ev, save=save_to_db)
    else:
        # quick orid/strike/dip/rake line
        print hp.output()
        p = 0
        
        if args.dbout:
            db = hp.output(format="ANTELOPE", dbout=args.dbout)
    # Done, return HashPype and/or FocalMechPlotter for debugging
    return hp, p


if __name__ == '__main__':
    hp, p = main()    

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

