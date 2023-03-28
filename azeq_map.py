from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
from obspy import read, Stream
from datetime import date
import matplotlib as mpl
mpl.use('macosx')

def _fig0(lat_ref, lon_ref, lats, lons, path, grid=True, save=False):
    """
    Plot azimuthal equidistant map centered on reference coordinates
    :param lat_ref: reference latitude
    :param lon_ref: reference longitude
    :param lats: Numpy Array of latitudes
    :param lons: Numpy Array of longitudes
    :param grid: draw grid at specified intervals. Defaults to true
    :param save: save figures. Defaults to false
    :param path: filepath to save figure
    :return: azimuthal equidistant map centered on reference coordinates
    """
    m = Basemap(projection='aeqd', lat_0=lat_ref, lon_0=lon_ref)
    m.drawmapboundary(fill_color='aqua')
    # draw coasts and fill continents.
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color='coral', lake_color='aqua')

    # 20 degree graticule
    if grid is True:
        m.drawparallels(np.arange(-80, 81, 20))
        m.drawmeridians(np.arange(-180, 180, 20))
    else:
        m.drawparallels()
        m.drawmeridians()

    # draw a black dot at the center.
    xpt, ypt = m(lon_ref, lat_ref)
    m.plot([xpt], [ypt], 'ko')
    # draw a black dot at event locations.
    xpts, ypts = m(lons, lats)
    m.plot([xpts], [ypts], 'ko')
    # draw the title.
    plt.title('Azimuthal Equidistant Map of Event Locations')

    if save is True:
        name = date.today()
        plt.savefig(os.path.join(path, 'azeq_{}.pdf'.format(name)))

    plt.show()

def check_data(filepath):
    """

    :param filepath: filepath of coordinates
    :return: data exists, true. Does not exist, false
    """
    data = np.loadtxt(filepath)
    if len(data[0]) > 1:
        return True
    else:
        return False

def get_event_data(enter_dir):
    """

    :param enter_dir: filepath of SAC files to plot events
    :return: event coordinates
    """
    evlats = []
    evlons = []
    if os.path.isdir(enter_dir):
        for path in os.listdir(enter_dir):
            if str(path[-3:]) == "sac":
                try:
                    st = read(os.path.join(enter_dir, path))
                    evlat = str(st[0].stats.sac.evla)
                    evlon = str(st[0].stats.sac.evlo)
                    evlats = np.append(evlats, evlat)
                    evlons = np.append(evlons, evlon)
                    # suppress errors so code doesn't stop
                except FileNotFoundError:
                    print("error for file: ", path)
            else:
                print("not a sac file")

    filename = input("Enter the name of the coordinate file:")
    figname = str(filename) + ".txt"
    with open("{}".format(figname), 'w+') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(evlats, evlons))

    return evlats, evlons, filename


# 0 if you don't have a file
filepath = '/Users/madeleinetan/Research/arrays_downloaded_112422/TA.C40A/test3.txt'

# column number (0 - N columns), 100 if you don't have one
lat_loc = 0

# column number (0 - N columns), 100 if you don't have one
lon_loc = 1
ref_lat = 42.28
ref_lon = 83.78

if filepath == 0:
    enter_dir = input("Enter the path of events: ")
    evlats, evlons, filename = get_event_data(enter_dir)
    print("Making map...")
    _fig0(ref_lat, ref_lon, evlats, evlons, enter_dir, grid=True, save=True)
else:
    loc = check_data(filepath)
    if loc is True:
        print("Data is in proper format")
        enter_dir = input("Enter the path of events: ")
        print("Making map...")
        data = np.loadtxt(filepath)
        _fig0(ref_lat, ref_lon, data[:,lat_loc], data[:,lon_loc], enter_dir, grid=True, save=True)
    else:
        print("Could not make map, data is only 1 column")
