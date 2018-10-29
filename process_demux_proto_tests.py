import general_tools
import gbw
import dumps
import energy_resol
import os

# ---------------------------------------------------------------------------
def process_demux_proto_tests(dirname):
    config = general_tools.get_conf()

    fulldirname = os.path.join(config['path_tests'], dirname)

    # -----------------------------------------------------------------------
    # Processing of hk files 
    hk = general_tools.get_hk(fulldirname, config)
    general_tools.plot_hk(hk, fulldirname, config)
 
    # -----------------------------------------------------------------------
    # Processing "BIAS, FEEDBAC and INPUT" dump files 
    dumps.process_dump(fulldirname, config, Max_duration=1.0)

    # -----------------------------------------------------------------------
    # Processing "Gain bandwidth characterization" 
    chan=0
    dumptype = "IQ-ALL"
    gbw.process_GBW(fulldirname, config, dumptype, chan)
    dumptype = "IQ-TST"
    gbw.process_GBW(fulldirname, config, dumptype, chan)

    # -----------------------------------------------------------------------
    # Processing "Carriers spectra characterization"
    pix=40 # test pixel index
    dumps.processIQ_multi(fulldirname, config, pix_zoom=pix, SR=False)

    # -----------------------------------------------------------------------
    # Processing "Energy resolution characterization"
    energy_resol.meas_energy_r(fulldirname, config)

# ---------------------------------------------------------------------------
