from DEMUX_TBox.general_tools import get_conf, print_conf, get_hk, print_hk, plot_hk
from DEMUX_TBox.gbw import process_GBW
from DEMUX_TBox.dumps import process_dump, processIQ_multi
from DEMUX_TBox.energy_resol import meas_energy_r
import os

# ---------------------------------------------------------------------------
def process_demux_proto_tests(dirname):
    config=get_conf()

    fulldirname = os.path.join(config['path_tests'], dirname)

    # -----------------------------------------------------------------------
    # Processing of hk files 
    hk=get_hk(fulldirname, config)
    plot_hk(hk, fulldirname, config)
 
    # -----------------------------------------------------------------------
    # Processing "BIAS, FEEDBAC and INPUT" dump files 
    process_dump(fulldirname, config, Max_duration=1.0)

    # -----------------------------------------------------------------------
    # Processing "Gain bandwidth characterization" 
    chan=0
    dumptype = "IQ-ALL"
    process_GBW(fulldirname, config, dumptype, chan)
    dumptype = "IQ-TST"
    process_GBW(fulldirname, config, dumptype, chan)

    # -----------------------------------------------------------------------
    # Processing "Carriers spectra characterization"
    pix=40 # test pixel index
    spt0dB = processIQ_multi(fulldirname, config, pix_zoom=pix, SR=False)

    # -----------------------------------------------------------------------
    # Processing "Energy resolution characterization"
    meas_energy_r(fulldirname, config)

# ---------------------------------------------------------------------------
