from obspy import Stream, read
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('macosx')

def _fig1(name, data_all, data45, data56, data67, data78):
    """
    Plot 5 stack traces on figure by station: all slowness values in s/km, [0.04-0.05), [0.05-0.06), [0.06-0.07), [0.07-0.08]
    Currently, data only accepts ObsPy format {STACKNAME}[0].data
    
    Parameters:
        name: title of figure; defaults to last 4 characters of subdirectory name
        data_all: stack of all data (trace shifted -0.15 Hz). Only accepts ObsPy format {STACKNAME}[0].data
        data45: stack of data [0.04-0.05) s/km. Only accepts ObsPy format {STACKNAME}[0].data
        data56: stack of data [0.05-0.06) s/km. Only accepts ObsPy format {STACKNAME}[0].data
        data67: stack of data [0.06-0.07) s/km. Only accepts ObsPy format {STACKNAME}[0].data
        data78: stack of data [0.07-0.08] s/km. Only accepts ObsPy format {STACKNAME}[0].data
        
    """
    plt.title('Stacking by slowness')
    plt.xlim([-5, 30])
    plt.ylim([-0.25, 0.6])
    plt.yticks([])
    # plt.grid(visible=True, which='major', axis='y', linewidth=0.5)
    plt.plot(time, data_all - 0.15, color='black', alpha=0.7)
    plt.plot(time, data45, color='blue', label='0.04-0.05')
    plt.plot(time, data56, color='green', label='0.05-0.06')
    plt.plot(time, data67, color='red', label='0.06-0.07')
    plt.plot(time, data78, color='orange', label='0.07-0.08')
    plt.legend(loc='best')
    title = str(name[-4:])
    plt.axhline(-0.15, -5, 30, color='black', alpha=0.5, lw=0.8)
    plt.axhline(0, -5, 30, color='black', alpha=0.5, lw=0.8)
    plt.savefig(os.path.join('/Users/madeleinetan/Research/stack_slow', '{}.png'.format(title)))
    plt.show(block=False)
    plt.pause(3)
    plt.close()

def deg2km(deg):
    """
    source: SeisPy
    
    Converts degrees to km
    
    Parameters:
        Deg: degrees
        
    """
    radius = 6371
    circum = 2*np.pi*radius
    conv = circum / 360
    km = deg * conv
    return km

def srad2skm(srad):
    """
    source: SeisPy
    
    Converts s/radian to s/km
    
    Parameters:
        srad: s/radian
        
    """
    sdeg = srad * ((2*np.pi)/360)
    return sdeg / deg2km(1)

# Main code ---
rootdir = '/Users/madeleinetan/Research/arrays_downloaded_112422' # data filepath
for j in os.listdir(rootdir): # j will be name input for _fig1
    enter_dir = os.path.join(rootdir, j)
    if os.path.isdir(enter_dir):
        # initialize streams
        st45 = Stream()
        st56 = Stream()
        st67 = Stream()
        st78 = Stream()
        stall = Stream()
        st = Stream()
        for f in glob.iglob(os.path.join(enter_dir, '*.sac')):
            st = read(f)
            stall += read(f)
            try:
                rayp = st[0].stats.sac.user4 # SAC ray parameter default container
            except AttributeError: # user4 may not be filled, if so calculates ray parameter manually
                evdp = st[0].stats.sac.evdp
                dist = st[0].stats.sac.dist
                arrivals = model.get_ray_paths(evdp, dist)
                arrival = arrivals[0]
                rayp = srad2skm(arrival.ray_param)

            # add traces by slowness to streams
            if 0.04 <= rayp < 0.05:
                st45 += read(f)
            elif 0.05 <= rayp < 0.06:
                st56 += read(f)
            elif 0.06 <= rayp < 0.07:
                st67 += read(f)
            elif 0.07 <= rayp < 0.08:
                st78 += read(f)
            else: # currently, traces outside [0.04-0.08] are not stacked separately (but are included in all stack)
                print('outside current bounds')

        # stack files within streams
        stack_all = stall.stack()
        stack45 = st45.stack()
        stack56 = st56.stack()
        stack67 = st67.stack()
        stack78 = st78.stack()

        # uses last SAC file in loop to find values
        tmin = st[0].stats.sac.b
        npts = st[0].stats.sac.npts
        dt = st[0].stats.sac.delta
        time = np.linspace(0, npts - 1, npts) * dt + tmin

        # set file canvas size
        f = plt.figure()
        f.set_figwidth(7)
        f.set_figheight(2)
        
        # plots traces and saves file
        _fig1(j, stack_all[0].data, stack45[0].data, stack56[0].data, stack67[0].data, stack78[0].data)
