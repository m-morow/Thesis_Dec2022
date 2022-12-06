#-- MT adapted based on Mijian Xu's Python code
# /*
#  * This code was taken from Mijian Xu's SeisPy
#  * SeisPy hk.py
#  *
#  * 2021
#  *
#  *
#  */

from obspy import read, Stream
from obspy.taup import TauPyModel
import glob
import numpy as np
import os
import csv
import matplotlib.pyplot as plt

def transarray(array, axis=0):
    if not isinstance(array, np.ndarray):
        raise ValueError('array should be `numpy.ndarray`')
    if len(array.shape) != 1:
        raise ValueError('array should be 1-d array')
    if axis == 0:
        return array.reshape(-1, array.shape[0])
    elif axis == 1:
        return array.reshape(array.shape[0], -1)
    else:
        raise ValueError('axis should be 0 or 1')

def vslow(v, rayp):
    return np.sqrt(1/(v**2) - rayp**2)


def tps(depth, eta_p, eta_s):
    return np.dot(transarray(eta_s - eta_p, axis=1), transarray(depth, axis=0))


def tppps(depth, eta_p, eta_s):
    return np.dot(transarray(eta_s + eta_p, axis=1), transarray(depth, axis=0))


def tpsps(depth, eta_s):
    return np.dot(transarray(2 * eta_s, axis=1), transarray(depth, axis=0))


def time2idx(times, ti0, dt):
    ti = ti0 + np.around(times / dt)
    return ti.reshape(ti.size).astype(int)

def hkstack(seis, t0, dt, p, h, kappa, vp=6.3, weight=(0.7, 0.2, 0.1)):
    # get dimensions
    nh = len(h)
    nk = len(kappa)
    nrf = len(p)

    # check the orientation of the seis array
    if seis.shape[0] != nrf:
        seis = seis.T
        if seis.shape[0] != nrf:
            raise IndexError('SEIS array dimensions should be (nt x nrf)')

    # amp correction for Ps
    am_cor = 151.5478 * p ** 2 + 3.2896 * p + 0.2618

    # get all vs, single column
    vs = vp / kappa

    # get index of direct P
    ti0 = round(t0 / dt)

    # initialize stacks
    tstack = np.zeros((nk, nh, 3))
    stack = np.zeros((nk, nh, 3))
    stack2 = np.zeros((nk, nh, 3))

    allstack = np.zeros((nk, nh, nrf))

    for i in range(nrf):
        eta_p = vslow(vp, p[i])
        eta_s = vslow(vs, p[i])

        # get times of Ps for all combinations of vs and H
        t1 = time2idx(tps(h, eta_p, eta_s), ti0, dt)
        t2 = time2idx(tppps(h, eta_p, eta_s), ti0, dt)
        t3 = time2idx(tpsps(h, eta_s), ti0, dt)

        tstack[:, :, 0] = am_cor[i] * seis[i, t1].reshape(nk, nh)
        tstack[:, :, 1] = am_cor[i] * seis[i, t2].reshape(nk, nh)
        tstack[:, :, 2] = -am_cor[i] * seis[i, t3].reshape(nk, nh)

        stack += tstack
        stack2 += tstack ** 2

        allstack[:, :, i] = weight[0] * tstack[:, :, 0] + weight[1] * tstack[:, :, 1] + weight[2] * tstack[:, :, 2]

    stack = stack / nrf
    stackvar = (stack2 - stack ** 2) / (nrf ** 2)

    allstackvar = np.var(allstack, axis=2)
    allstack = np.mean(allstack, axis=2)
    Normed_stack = allstack - np.min(allstack)
    Normed_stack = Normed_stack / np.max(Normed_stack)
    return stack, stackvar, Normed_stack, allstackvar, allstack

def _fig0(h, k, stack, besth, bestk):
    """
    Plot Hk stack
    :param h: array
    :param k: array
    :param stack: stacked Hk matrix (best to use Normed stack)
    :param besth: grab from ci output
    :param bestk: grab from ci output
    :return:
    """
    max_level = 1
    min_level = 0
    step_level = 0.005
    plt.title("H-k stack for all stations")
    plt.xlabel("H (km)")
    plt.ylabel("k")
    plt.xticks(ticks=[30, 35, 40, 45, 50, 55], labels=["30", "35", "40", "45", "50", "55"])
    plt.contourf(h, k, stack, levels=np.arange(min_level, max_level + step_level, step_level),
                 vmin=0, vmax=1, cmap='jet')
    lab = str(besth) + " km, " + str("{:.3f}".format(bestk))
    plt.scatter(besth, bestk, c='k', marker='+', s=70, label=lab)
    plt.legend(loc="best")

def ci(allstack, h, kappa, ev_num):
    """
    Search best H and kappa from stacked matrix.
    Calculate error for H and kappa
    :param allstack: stacked HK matrix
    :param h: 1-D array of H
    :param kappa: 1-D array of kappa
    :param ev_num: event number
    :return:
    """
    [i, j] = np.unravel_index(allstack.argmax(), allstack.shape)
    bestk = kappa[i]
    besth = h[j]

    cvalue = 1 - np.std(allstack.reshape(allstack.size)) / np.sqrt(ev_num)
    # cs = plt.contour(h, kappa, allstack, cvalue)
    # cs_path = cs.collections[0].get_paths()[0].vertices
    # maxhsig = (np.max(cs_path[:, 0]) - np.min(cs_path[:, 0])) / 2
    # maxksig = (np.max(cs_path[:, 1]) - np.min(cs_path[:, 1])) / 2
    plt.close()
    return besth, bestk, cvalue

# main code -------------------

rootdir = '/Users/madeleinetan/Research/TA_ARRAYS/T1/'
for j in os.listdir(rootdir):
    enter_dir = os.path.join(rootdir, j)
    if os.path.isdir(enter_dir):

        # initialize arrays
        st = Stream()
        suffix = str('*.sac')
        p = []

        ev_num = 1
        for f in glob.iglob(os.path.join(enter_dir, suffix)):
            # read in traces in filepath to stream
            st += read(f)
            ev_num += 1

        array = np.array(st)
        for i, tr in enumerate(st):
            evdp = st[i].stats.sac.evdp
            deg = st[i].stats.sac.gcarc
            try:
                p = np.append(p, st[i].stats.sac.user4)
            except AttributeError:
                model = TauPyModel(model="ak135")
                arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                                  distance_in_degree=deg,
                                                  phase_list=["P"])
                p = np.append(p, ((arrivals[0].ray_param * (1 / (180 * np.pi))) / 111.195))

        hmax = 55.5
        hmin = 30
        hstep = 0.5
        h = np.arange(hmin, hmax, hstep)
        kmax = 2.0
        kmin = 1.55
        kstep = 0.01
        k = np.arange(kmin, kmax, kstep)

        stack, stackvar, Normed_stack, allstackvar, allstack = hkstack(array, 10, 0.1, p=p, h=h,
                                                                       kappa=k, vp=6.3, weight=(0.7, 0.2, 0.1))

        besth, bestk, cvalue = ci(allstack, h=h, kappa=k, ev_num=ev_num)

        # plot figure
        _fig0(h, k, Normed_stack, besth, bestk)

        # save figure
        # figname = str(j)+".hk.pdf"
        figname = "T1.hk.pdf"
        plt.savefig(os.path.join('/Users/madeleinetan/Research/TA_ARRAYS/figures_112022/hk/', figname))
        plt.close()
