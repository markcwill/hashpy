#
import os

from obspy.core.event import read_events
from obspy.io.quakeml.core import Pickler, Catalog

from hashpy.hashpype import HashPype
import hashpy.io.obspy as opio

PWD = os.path.dirname(__file__)

QMLFILE = PWD+"/data/nn00450180_phase.xml"
VZFILE = PWD+"/data/vz.pickema2"


def test_rid():
    opio.PREFIX = "smi:local.test-server"
    ri = opio.rid("RID")
    assert str(ri) == "smi:local.test-server/RID"


def test_obspyinput():
    cat = read_events(QMLFILE)
    assert len(cat.events) > 0
    hp = HashPype(
        input_factory=opio.input_event,
        output_factory=opio.output_event,
    )
    hp.input(cat.events[0])
    assert hp.npol == 34


def test_full_integration():
    """Full integration test"""
    #
    # TODO: hard code in HASH default settings to this test in case they change
    #
    cat = read_events(QMLFILE)
    assert len(cat.events) > 0
    hp = HashPype(
        input_factory=opio.input_event,
        output_factory=opio.output_event,
    )
    hp.input(cat.events[0])
    hp.vmodels = [VZFILE]
    hp.driver2(False, False) # this breaks py.test
    fmevent = hp.output()
    
    # assert equality of expected solution parts: picks, maybe str/dip/rake if
    # deterministic, but HASH is weird that way.
    assert len(fmevent.picks) == 34
    assert hp.qual[0] == "B"
    qml = Pickler().dumps(Catalog(events=[fmevent]))
    print qml

