import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
def read_events(eventsfilename):
       
    dt2=np.dtype([('timestamp', np.float), \
                 ('channelId', np.int8), \
                 ('pixelId', np.int8), \
                 ('energy', np.float32), \
                 ('offset', np.float32)])
                                  
    fdat=open(eventsfilename, 'rb')
    event_list=np.fromfile(fdat, dtype=dt2)
    fdat.close()
    return(event_list[:]['timestamp'], event_list[:]['channelId'], event_list[:]['pixelId'], event_list[:]['energy'], event_list[:]['offset'])
# -----------------------------------------------------------------------------

fname='20190520-112338_check_resol.events'
t,chid,pixid,e,b=read_events(fname)

print(e.min(), e.max())
print(len(e))

plt.plot(t,e)
plt.show()