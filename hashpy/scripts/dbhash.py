#
"""
hashpy.scripts.dbhash

- Mark Williams (2013) 
- Nevada Seismological Laboratory

NSL's production CLI HASH program: BEWARE!!

Example program to run HASH using antelope
    This is esoteric to NSL, but pretty straighforward to modify. For
    example, this version is hard-coded to not enforce minimum gap check.

    In graphics mode, this program may create a folder and/or save an
    image file. 

Call by importing main() into external executable script
"""
import sys
import os
import logging
from argparse import ArgumentParser
from hashpy.hashpype import HashPype, HashError
from hashpy.io.antelopeIO import load_pf, eventfocalmech2db, dbloc_source_db

LOG = logging.getLogger()

parser = ArgumentParser()
parser.add_argument("dbin", help="Input database")
parser.add_argument("dbout", help="Output database", nargs='?')
parser.add_argument("-p", "--plot", help="Plot result", action='store_true')
parser.add_argument("-l", "--loc", help="dbloc2 mode", action='store_true')
parser.add_argument("-i", "--image", help="Save image with db", action='store_true')
parser.add_argument("--pf",   help="Parameter file")
group = parser.add_mutually_exclusive_group() #required=True)
group.add_argument("--evid", help="Event ID", type=int)
group.add_argument("--orid", help="Origin ID", type=int)


def run_hash(dbname, orid=None, pf=None):
    """
    Perform a 'driver2' HASH run using database input.

    Input
    -----
    dbname : str of db name
    orid : int of ORID
    pf : str of pf file name
    
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
    
    # Run driver script method
    hp.driver2(check_for_maximum_gap_size=False)
    return hp
    

class SaveFunction(object):
    """
    Create a save function for matplotlib GTK Figure
    """
    @staticmethod
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
    
    def __init__(self, dbname, dump_bitmap):
        self.dbname = dbname
        self.dump_bitmap = dump_bitmap

    def __call__(self, fmplotter):
        focal_mech = fmplotter.event.focal_mechanisms[fmplotter._fm_index]
        if focal_mech is not fmplotter.event.preferred_focal_mechanism():
            fmplotter.event.preferred_focal_mechanism_id = str(focal_mech.resource_id)
        # Save to db        
        eventfocalmech2db(event=fmplotter.event, database=self.dbname)
        if self.dump_bitmap:
            vers = fmplotter.event.preferred_origin().creation_info.version
            dbdir = os.path.dirname(self.dbname)
            self._dump_bitmap(figure=fmplotter.fig, directory=dbdir, uid=vers)


def dbhash(args):
    """
    Run HASH using database given a namespace of arguments
    """
    # Special 'dbloc2' settings
    if args.loc:
        import curds2.dbapi2 as dbapi2
        from curds2.cursors import InteractiveCursor
        # alter args b/c dbloc2 passes a db and a row number
        args.dbin = args.dbin.rstrip('.origin')
        curs = dbapi2.connect(args.dbin, cursor_factory=InteractiveCursor).cursor()
        n = curs.execute.lookup(table='origin')
        rec = int(args.dbout)
        curs.scroll(rec, 'absolute')
        args.orid = curs.fetchone()['orid']
        args.dbout = dbloc_source_db(args.dbin, pointer=False)
        args.plot = True   # force plot
        args.image = True  # force saving image to db folder
        curs.close()

    #--- Run HASH ----------------------------------------------------#
    hp = run_hash(args.dbin, orid=args.orid, pf=args.pf)

    #--- OUTPUT --- Launch plotter or spit out solution --------------#
    if args.plot:
        from hashpy.plotting.focalmechplotter import FocalMechPlotter
        save_plot_to_db = SaveFunction(args.dbout, args.image)
        ev = hp.output(format="OBSPY")
        p = FocalMechPlotter(ev, save=save_plot_to_db)
    else:
        # quick orid/strike/dip/rake line
        print(hp.output())
        p = 0    
        if args.dbout:
            db = hp.output(format="ANTELOPE", dbout=args.dbout)
    # Done, return HashPype and/or FocalMechPlotter for debugging
    return hp, p


def main():
    """
    CLI program to run dbhash
    
    Handle input args, call function, handle output.
    """
    # Root logger to stderr for now
    logging.basicConfig(format='[%(levelname)s]: %(message)s')
    
    try:
        args = parser.parse_args()
        out = dbhash(args)
        return 0
    except Exception as e:
        LOG.exception(e)
        LOG.critical("Uncaught error, exiting dbhash")
        return 1


#-- Testing only -----------------------------------------------------#
if __name__ == '__main__':
    ret = main()
