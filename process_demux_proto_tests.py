import general_tools
import hk_tools
import gbw
import dumps
import dumps_dds
import baseline
import os
import ep_tools
import scan_feedback_tools

# ---------------------------------------------------------------------------
def process_demux_proto_tests(dirname, verbose=False):

    # -----------------------------------------------------------------------
    # Reading demux and session informations 
    config = general_tools.get_csv('demux_tools_cfg.csv')

    fulldirname = os.path.join(os.path.normcase(config['path_tests']), dirname)
    session_info = general_tools.get_csv(os.path.join(fulldirname, config['session_info']))

    if verbose:
        general_tools.print_dict(session_info, 'demux')
        general_tools.print_dict(session_info, 'session')

    # -----------------------------------------------------------------------
    # Processing of hk files 
    hk, hk_lims = hk_tools.get_hk(fulldirname, config)
    if hk != 0:
        hk_tools.plot_hk(hk, hk_lims, fulldirname, config, plt_temp=True)
 
    # -----------------------------------------------------------------------
    # Processing scan feedback data 
    sfb_dat = scan_feedback_tools.get_scanfb(fulldirname, config)
    if sfb_dat !=0:
        scan_feedback_tools.plot_scanfb(sfb_dat, fulldirname, config)

    # -----------------------------------------------------------------------
    # Processing "Carriers spectra characterization"
    tst_pix=40 # test pixel index
    _, pix_pos=dumps.process_iq_multi(fulldirname, config, pix_zoom=tst_pix)
    dumps.process_iq_tst_multi(fulldirname, config, window=False, bw_correction=True)

    # -----------------------------------------------------------------------
    # Processing "BIAS, FEEDBAC and INPUT" dump files 
    dumps.process_dump(fulldirname, config, max_duration=1.0, pix_id=pix_pos)
    dumps_dds.process_dump_dds(fulldirname, config, max_duration=1.0)

    # -----------------------------------------------------------------------
    # Processing "Gain bandwidth characterization" 
    channel=0
    gbw.process_gbw(fulldirname, config, channel)

    # -----------------------------------------------------------------------
    # Checking baseline level
    baseline.check_baseline(fulldirname, config)

    # -----------------------------------------------------------------------
    # Checking delock behaviour
    dumps.process_dump_nl(fulldirname, config)
    dumps.process_dump_delock_iq(fulldirname, config, "NL_anti-Delock-OFF")
    dumps.process_dump_delock_iq(fulldirname, config, "NL_anti-Delock--ON")

    # -----------------------------------------------------------------------
    # Checking pulse generator behaviour
    dumps.process_dump_pulses_adc_dac(fulldirname, config, 'IN-BIA_PULSE', zoom_factor=50)
    dumps.process_dump_pulses_adc_dac(fulldirname, config, 'IN-FBK_PULSE', zoom_factor=50)
    dumps.process_dump_pulses_iq(fulldirname, config)

    # -----------------------------------------------------------------------
    # Measuring energy resolution
    ep_tools.ep(fulldirname, config)

# ---------------------------------------------------------------------------
