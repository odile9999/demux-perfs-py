import general_tools
import hk_tools
import gbw
import dumps
import energy_resol
import baseline
import os

# ---------------------------------------------------------------------------
def process_demux_proto_tests(dirname):
    config = general_tools.get_conf()

    fulldirname = os.path.join(os.path.normcase(config['path_tests']), dirname)

    # -----------------------------------------------------------------------
    # Processing of hk files 
    hk, hk_lims = hk_tools.get_hk(fulldirname, config)
    if hk != 0:
        hk_tools.plot_hk(hk, hk_lims, fulldirname, config)
 
    # -----------------------------------------------------------------------
    # Processing "BIAS, FEEDBAC and INPUT" dump files 
    #dumps.process_dump(fulldirname, config, Max_duration=1.0)

    # -----------------------------------------------------------------------
    # Processing "Gain bandwidth characterization" 
    channel=0
    dumptype = "IQ-ALL"
    gbw.process_GBW(fulldirname, config, dumptype, channel)
    dumptype = "IQ-TST"
    gbw.process_GBW(fulldirname, config, dumptype, channel)

    # -----------------------------------------------------------------------
    # Processing "Carriers spectra characterization"
    pix=40 # test pixel index
    dumps.processIQ_multi(fulldirname, config, pix_zoom=pix)
    dumps.processIQ_TST_multi(fulldirname, config, window=False, BW_CORRECTION=True)

    # -----------------------------------------------------------------------
    # Checking baseline level
    baseline.check_baseline(fulldirname, config)

    # -----------------------------------------------------------------------
    # Processing "Energy resolution characterization"
    energy_resol.meas_energy_r(fulldirname, config)

# ---------------------------------------------------------------------------
