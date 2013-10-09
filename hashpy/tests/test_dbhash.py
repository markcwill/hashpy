#!/usr/bin/env python
"""
#  dbhash.py
#
# functions which run HASH using Antelope and hashpy
#
# optionally, can plot by passing ObsPy output to a plotter
"""
import os
from hashpy.hashpype import HashPype, HashError
from hashpy.io.antelopeIO import ( load_pf, eventfocalmech2db, dbloc_source_db, RowPointerDict )

def _dump_bitmap(figure=None, directory='.', uid=None):
    """
    Dump mpl figure to a file
    
    Given a figure & dir, make folder called 'images' and
    dump it into a png file called '<uid>_focalmech.png'

    """
    filename = 'focalmech.png'
    if uid:
        filename = '_'.join([str(uid),filename])
    imagedir = os.path.join(directory, 'images')
    if not os.path.exists(imagedir):
        os.mkdir(imagedir)
    fullname = os.path.join(imagedir, filename)
    figure.savefig(fullname)

def dbhash_run(dbname, orid=None, pf=None):
    """
    Perform a HASH run using Database input and command line args
    
    Input
    -----
    args : Namespace of command line args
    
    Returns : hp : a hashpy.HashPype object containing solutions

    """
    hp = HashPype()
    
    # Load settings data from a pf file...
    if pf:
        load_pf(hp, pffile=pf)
    else:
        load_pf(hp)
    
    # Grab data from the db...
    hp.input(dbname, format="ANTELOPE", orid=orid)
    
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
    parser.add_argument("-i", "--image", help="Save image with db", action='store_true')
    parser.add_argument("--pf",   help="Parameter file")
    group = parser.add_mutually_exclusive_group() #required=True)
    group.add_argument("--evid", help="Event ID", type=int)
    group.add_argument("--orid", help="Origin ID", type=int)
    args = parser.parse_args()
    
    # Special 'dbloc2' settings
    if args.loc:
        from antelope.datascope import Dbptr
        # alter args b/c dbloc2 passes a db and a row number
        args.dbin = args.dbin.rstrip('.origin')
        db = Dbptr(args.dbin)
        db = db.lookup(table='origin')
        db.record = int(args.dbout)
        args.orid = db.getv('orid')[0]
        args.dbout = dbloc_source_db(args.dbin, pointer=False)
        args.plot = True   # force plot
        args.image = True  # force saving image to db folder
    
    # Now that we have a save location from command line args,
    # make a function to save to that database. The plotter is I/O
    # agnostic, it will accept a function to save anything anyhow anywhichway
    #
    def save_plot_to_db(fmplotter, dbname=args.dbout, dump_bitmap=args.image):
        focal_mech = fmplotter.event.focal_mechanisms[fmplotter._fm_index]
        if focal_mech is not fmplotter.event.preferred_focal_mechanism():
            fmplotter.event.preferred_focal_mechanism_id = focal_mech.resource_id.resource_id
        
        # Save to db        
        eventfocalmech2db(event=fmplotter.event, database=dbname)
        
        if dump_bitmap:
            vers = fmplotter.event.preferred_origin().creation_info.version
            dbdir = os.path.dirname(dbname)
            _dump_bitmap(figure=fmplotter.fig, directory=dbdir, uid=vers)

    # Run HASH
    hp = dbhash_run(args.dbin, orid=args.orid, pf=args.pf)

    # Launch plotter or spit out solution
    if args.plot:
        from hashpy.plotting.focalmechplotter import FocalMechPlotter
        ev = hp.output(format="OBSPY")
        p = FocalMechPlotter(ev, save=save_plot_to_db)
    else:
        # quick orid/strike/dip/rake line
        print hp.output()
        p = 0    
        if args.dbout:
            db = hp.output(format="ANTELOPE", dbout=args.dbout)
    
    # Done, return HashPype and/or FocalMechPlotter for debugging
    return hp, p


if __name__ == '__main__':
    import time
    t0 = time.time()
    ret = main()
    print time.time() - t0
