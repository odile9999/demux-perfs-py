#
# This is a python script to be used to process the data of specific test sequences.
# 

from process_demux_proto_tests import process_demux_proto_tests

#dirname = '20181017_091738'
#dirname = '20181019_160619'
#dirname = '20181115_095600_fw0000'
#dirname = '20181115_101500_fwv16'
#dirname = '20181115_103300_fw16opt'
#dirname = '20181115_124600_fwv19'
#dirname = '2018-11-15_17.53.30__perfo-proto'
#dirname = '2018-11-16_09.32.37__perfo-proto'
#dirname = '2018-11-16_09.48.39__perfo-proto'
#dirname = '2018-11-16_10.04.10__perfo-proto'
#dirname = '2018-11-16_11.52.34__no_spread'
#dirname = '2018-11-16_16.44.40__FB_trunc=2'
#dirname = '2018-11-16_17.08.47__FB_trunc=4'
dirname = '2018-11-19_16.41.35__validation du fw0100'

process_demux_proto_tests(dirname)
