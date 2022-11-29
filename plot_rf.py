import matplotlib.pyplot as plt
from matplotlib import gridspec
from obspy import read, Stream
import numpy as np
import os
import glob
import matplotlib as mpl
mpl.use('macosx')


def _fig0(data, scale_factor=1, offset, time):
    """
    Plot receiver functions as a function of back-azimuth (deg)

    Parameters:
        data: np.array
        scale_factor: scale trace up (f > 1), down (0 < f < 1). Default f=1
        offset: space traces out
        time: y-axis, calculated from sac header

    """
    ax0.set_ylim([25, 0])
    ax0.set_xlim([0, 380])
    ax0.set_xticks([])
    ax0.set_ylabel('time lag, s')
    ax0.minorticks_on()
    ax0.grid(which='major', axis='y', linewidth=0.5)
    fig.suptitle('TA.{}_M6.5+'.format(station))
    ax0.plot(data * scale_factor + offset, time, color='black', linewidth=0.8)
    ax0.fill_betweenx(time, data * scale_factor + offset, offset, where=(data * scale_factor + offset) > offset,
                      color='r', alpha=0.7)
    ax0.fill_betweenx(time, data * scale_factor + offset, offset, where=(data * scale_factor + offset) < offset,
                      color='b', alpha=0.7)

def _fig1(offset, backazimuth):
    """
    Scatter plot of receiver function distances

    Parameters:
        offset: space points out
        backazimuth: back-azimuth (deg)

    """
    ax1.set_xlim(0, 380)
    ax1.set_ylim(0, 360)
    ax1.set_ylabel('back azimuth')
    ax1.set_xticks([])
    ax1.scatter(offset, backazimuth, c='k', s=10)

def _fig2(data, time, scale_factor=1):
    """
    Plots receiver function stack [0, 25] seconds

    Parameters:
        data: np.array
        time: calculated using sac header
        scale_factor: scale trace up (f > 1), down (0 < f < 1). Default f=1

    """
    ax2.set_ylim([25, 0])
    ax2.grid(visible=True, which='major', axis='y', linewidth=0.5)
    ax2.set_yticks([])
    ax2.plot(data*scale_factor, time, color='black', linewidth=0.8)
    ax2.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) > 0, color='r', alpha=0.7)
    ax2.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) < 0, color='b', alpha=0.7)

def _fig3(data, time, scale_factor=1):
    """
    Plots receiver function stack [0, 70] seconds

    Parameters:
        data: np.array
        time: calculated using sac header
        scale_factor: scale trace up (f > 1), down (0 < f < 1). Default f=1

    """
    ax3.set_ylim([70, 0])
    ax3.yaxis.tick_right()
    ax3.plot(data*scale_factor, time, color='black', linewidth=0.8)
    ax3.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) > 0, color='r', alpha=0.7)
    ax3.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) < 0, color='b', alpha=0.7)

rootdir = '/Users/madeleinetan/Research/TA_ARRAYS/arrays_downloaded_112022/'
for j in os.listdir(rootdir):
    enter_dir = os.path.join(rootdir, j)
    if os.path.isdir(enter_dir):

        # initialize arrays
        backaz = []
        fff = []
        filenames = []
        st = Stream()
        suffix = str('*.sac')

        for f in glob.iglob(os.path.join(enter_dir, suffix)):
            # read in traces in filepath to stream
            st += read(f)
            filenames = np.append(filenames, f)

        for tr in st:
            # get back-azimuth of trace
            backaz = np.append(backaz, tr.stats.sac.baz)

        # sort by increasing back-azimuth
        idx = np.argsort(backaz)

        st1 = Stream()

        for file in filenames[idx]:
            # index filenames by increasing back-azimuth
            filepath = '{}'.format(file)
            for f in glob.iglob(filepath):
                # read trace into stream by increasing back-azimuth
                st1 += read(f)

        # set parameters using first trace in stream
        station = str(st[0].stats.station)
        tmin = st1[0].stats.sac.b
        npts = st1[0].stats.sac.npts
        dt = st1[0].stats.sac.delta
        time = np.linspace(0, npts - 1, npts) * dt + tmin

        # set parameters using first trace in stream
        fig = plt.figure()
        grid = plt.GridSpec(4, 5, hspace=0.2, wspace=0.2)
        ax0 = fig.add_subplot(grid[0:3, 0:3])  # rfs
        ax1 = fig.add_subplot(grid[3, 0:3])  # baz plot
        ax2 = fig.add_subplot(grid[0:3, 3])  # stack
        ax3 = fig.add_subplot(grid[0:4, 4])  # extended stack

        # calculate offset of traces in record section
        offset = 180 - (((len(idx)) / 2) * 10)
        plot_offset = 10
        if offset < 0:
            print("WARNING: OFFSET IS NEGATIVE!")
            plot_offset = 5


        for j, tr in enumerate(st1):
            backazimuth = tr.stats.sac.baz

            # plot RF trace
            _fig0(tr.data, 30, offset, time)

            # scatter plot of trace back-azimuth
            _fig1(offset, backazimuth)

            offset = offset + plot_offset

        # stack traces in stream st1
        stack = st1.stack()

        for k, trr in enumerate(stack):

            # plot RF stack [0, 25] s
            _fig2(trr.data, time)

            # plot RF [0, 70] s
            _fig3(trr.data, time)

        # save figures as pdf
        plt.savefig('/Users/madeleinetan/Research/TA_ARRAYS/figures_112022/{}.pdf'.format(station))

        # close and clear figure to decrease memory space
        plt.close()
        print("success for:", enter_dir)
