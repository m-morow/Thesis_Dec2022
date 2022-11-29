import matplotlib.pyplot as plt
from matplotlib import gridspec
from obspy import read, Stream, Trace
import numpy as np
import os
import glob
import matplotlib as mpl
mpl.use('macosx')

def _fig0(data, scale_factor=1, offset, time):
    """
    Plot receiver functions as a function of distance (deg)

    Parameters:
        data: np.array
        scale_factor: scale trace up (f > 1), down (0 < f < 1). Default f=1
        offset: space traces out
        time: y-axis, calculated from sac header

    """
    ax0.set_ylim([80, 0])
    ax0.set_xlim([10, 110])
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

def _fig1(offset, distance):
    """
    Scatter plot of receiver function distances

    Parameters:
        offset: space points out (follows distance values)
        distance: distance (deg)

    """
    ax1.set_xlim(10, 110)
    ax1.set_ylim(10, 110)
    ax1.set_ylabel('distance')
    ax1.set_xticks([])
    ax1.scatter(offset, distance, c='k', s=10)

def _fig2(data, time, scale_factor=1):
    """
    Plots receiver function stack

    Parameters:
        data: np.array
        time: calculated using sac header
        scale_factor: scale trace up (f > 1), down (0 < f < 1). Default f=1

    """
    ax2.set_ylim([80, 0])
    ax2.set_yticks([])
    ax2.plot(data*scale_factor, time, color='black', linewidth=0.8)
    ax2.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) > 0, color='r', alpha=0.7)
    ax2.fill_betweenx(time, data*scale_factor, x2=0, where=(data*scale_factor) < 0, color='b', alpha=0.7)

rootdir = '/Users/madeleinetan/Research/TA_ARRAYS/arrays_downloaded_112022/'
for j in os.listdir(rootdir):
    enter_dir = os.path.join(rootdir, j)
    if os.path.isdir(enter_dir):

        # initialize arrays
        dist = []
        distances = []
        delta_t = []
        st_moveout = []
        filenames = []
        st = Stream()

        suffix = str('*.sac')

        for f in glob.iglob(os.path.join(enter_dir, suffix)):
            # read in traces in filepath to stream
            st += read(f)
            filenames = np.append(filenames, f)

        for tr in st:
            # get distance of trace
            distances = np.append(distances, tr.stats.sac.dist)

        # determine average distance to use as reference distance
        ref_dist = np.average(distances / 111.195)
        # sort by increasing distance
        idx = np.argsort(distances)

        st1 = Stream()

        for file in filenames[idx]:
            # index filenames by increasing distance
            file_path = '{}'.format(file)
            for f in glob.iglob(file_path):
                # read trace into stream by increasing distance
                st1 += read(f)

        for tr in st1:
            dist = tr.stats.sac.dist
            # calculate t based on simple linear moveout
            # t [s] = slowness [s/deg] * (trace distance [deg] - reference distance)
            t = 0.07 * ((dist / 111.195) - ref_dist)
            deltat = np.append(delta_t, t)

        for i, tr in enumerate(st1):
            # convert obspy trace to np.array
            tr2 = np.array(tr.data)
            # determine how many shifts in delta t you need (interval = 0.1 s)
            shift = int(np.rint(delta_t[i]/0.1))
            # move elements by the shift value you determined
            tr3 = np.roll(tr2, shift)
            # convert to obspy trace
            tr4 = Trace(data=tr3)
            # read traces into st_moveout stream
            st_moveout += Stream(traces=[tr4])

        # set parameters using first trace in stream
        station = str(st[0].stats.station)
        tmin = st[0].stats.sac.b
        npts = st[0].stats.sac.npts
        dt = st[0].stats.sac.delta
        time = np.linspace(0, npts - 1, npts) * dt + tmin

        # initialize figure and subplots
        fig = plt.figure()
        grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
        ax0 = fig.add_subplot(grid[0:3, 0:3])  # rfs
        ax1 = fig.add_subplot(grid[3, 0:3])  # baz plot
        ax2 = fig.add_subplot(grid[0:3, 3])  # stack

        for j, tr in enumerate(st1):
            dist = tr.stats.sac.dist
            offset = dist/111.195

            # plot RF trace
            _fig0(tr.data, 15, offset, time)

            # scatter plot of trace distance
            _fig1(offset, offset)


        # plot RF stack
        stack = np.sum([tr.data for tr in st_moveout], axis=0)
        _fig2(stack, time, 1)

        # save figure as pdf
        plt.savefig('/Users/madeleinetan/Research/TA_ARRAYS/figures_112022/slant_stacks/{}.pdf'.format(station))

        # close and clear figure to decrease memory space
        plt.close()
        print("success for:", enter_dir)


"""Check delta t with data"""
# plt.title("\u0394 t vs. Distance")
# plt.xlabel("Distance from event (deg)")
# plt.ylabel("\u0394 t (s)")
# plt.scatter(dist/111.195, deltat, c='black')
# plt.scatter(58.375, 0, marker="*", s=50, c='red')
# plt.savefig('/Users/madeleinetan/Research/TA_ARRAYS/figures_112022/TA.C40A_deltat.pdf')
# plt.show()


# counts, bins = np.histogram(dist/111.195)
# plt.hist(bins[:-1], bins, weights=counts)
# plt.show()
# avg_dist = np.average(dist/111.195)
# print(avg_dist)

