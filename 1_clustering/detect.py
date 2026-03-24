import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator, LongitudeLocator)

import matplotlib.pyplot as plt
import shapefile
import numpy as np
import pandas as pd
from cmcrameri import cm
import obspy
import matplotlib.patheffects as path_effects
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84

import requests
from bs4 import BeautifulSoup
from netCDF4 import Dataset
import subprocess
import shlex
from matplotlib.colors import LightSource
import matplotlib as mpl
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN

import matplotlib.dates as mdates

from obspy.io.nied import knet
import obspy

#datarootdir = '/Users/onion/GoogleDrive/Work/2024Noto/submission_contents_2024Noto/github_repository/materials'
datarootdir = '/Users/onion/GoogleDrive/Work/2024Noto/submission_contents_2024Noto/2024Noto/materials'


def map_tickslabels(ax, xtickint, ytickint, lonmin, lonmax, latmin, latmax):
    xticklabels = np.arange(lonmin - lonmin%xtickint, lonmax - lonmax%xtickint + xtickint, xtickint)
    yticklabels = np.arange(latmin - latmin%ytickint, latmax - latmax%ytickint + ytickint, ytickint)
    ax.set_xticks(xticklabels, crs=ccrs.PlateCarree())
    ax.set_yticks(yticklabels, crs=ccrs.PlateCarree())
    if isinstance(xtickint, int):
        lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.0f')
    else:
        lon_formatter = LongitudeFormatter(zero_direction_label=True, number_format='.'+str(len(str(xtickint).replace('.',' ').split()[1]))+'f')
    if isinstance(ytickint, int):
        lat_formatter = LatitudeFormatter(number_format='.0f')
    else:
        lat_formatter = LatitudeFormatter(number_format='.'+str(len(str(ytickint).replace('.',' ').split()[1]))+'f')
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent([lonmin, lonmax, latmin, latmax])
    return ax

from cycler import cycler
#cmap2 = plt.get_cmap('Pastel1', 8)
cmap2 = plt.get_cmap('Set2', 8)
custom_color_cycle=[]
for i in range(cmap2.N):
    rgb = cmap2(i)[:3]
    custom_color_cycle.append(str(mpl.colors.rgb2hex(rgb)))
cyc_pastel = cycler(color=custom_color_cycle)
#ax.set_prop_cycle(cyc_pastel)
plt.rc('axes', prop_cycle=(cycler(color=custom_color_cycle)))

import matplotlib as mpl
initfontsize = 10
mpl.rc('axes', labelsize=initfontsize, titlesize=initfontsize)
mpl.rc('xtick', labelsize=initfontsize)
mpl.rc('ytick', labelsize=initfontsize)
mpl.rc('legend', fontsize=initfontsize, edgecolor='none')
mpl.rc('savefig', dpi=600, transparent=False)
mpl.rc('font', size=initfontsize)

mpl.rcParams['font.weight'] = 400
mpl.rcParams['font.family'] = 'Open Sans'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Open Sans'
mpl.rcParams['mathtext.it'] = 'Open Sans:italic'
mpl.rcParams['mathtext.bf'] = 'Open Sans:bold'

import zipfile
import os
import io
from netCDF4 import Dataset
from fastkml.kml import KML

red=plt.get_cmap('RdBu_r')(220)


import obspy
from pyrocko import moment_tensor as pmt
#from pyrocko.plot import beachball
from obspy.imaging.mopad_wrapper import beach


# load model parameters from `fort.40` and store the values as pandas data frame
def load_fort40(infile):
    col_names = [ 'c{0:02d}'.format(i) for i in range(10)]
    df = pd.read_table(infile, names=col_names, header=None, delimiter='\t')
    tmp1 = [ float(x) for x in df['c00'][1].replace(' ', ',').split(',') if x ]
    tmp3 = [ float(x) for x in df['c00'][3].replace(' ', ',').split(',') if x ]
    tmp5 = [ float(x) for x in df['c00'][5].replace(' ', ',').split(',') if x ]
    tmp7 = [ float(x) for x in df['c00'][7].replace(' ', ',').split(',') if x ]
    df = pd.DataFrame({ 'moment'   : tmp1[0],
                        'mw'       : tmp1[1],
                        'rigidity' : tmp1[2],
                        'lat'      : tmp1[3],
                        'lon'      : tmp1[4],
                        'depth'    : tmp1[5],
                        'vr'       : tmp1[6],
                        'nsurface' : tmp1[7],

                        'strike'   : tmp3[0],
                        'dip'      : tmp3[1],

                        'xx'       : tmp5[0],
                        'yy'       : tmp5[1],
                        'mn'       : int(tmp5[2]),
                        'nn'       : int(tmp5[3]),
                        'm0'       : int(tmp5[4]),
                        'n0'       : int(tmp5[5]),
                        'tr'       : tmp5[6],
                        'jtn'      : int(tmp5[7]),
                        'icmn'     : int(tmp5[8]),

                        'variance' : tmp7[0],

                      }, index=[0])
    return df

# select a preferred fault plane
def selectplane(modelstk, modeldip, stk0, dip0, rake0, stk1, dip1, rake1):
    vecmodelplane = faultnormalvec(modelstk, modeldip)
    vecplane0 = faultnormalvec(stk0, dip0)
    vecplane1 = faultnormalvec(stk1, dip1)
    tmp0 = np.inner(vecmodelplane, vecplane0)
    tmp1 = np.inner(vecmodelplane, vecplane1)
    if abs(tmp0) > abs(tmp1):
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    elif abs(tmp0) < abs(tmp1):
        stk_s = stk1
        dip_s = dip1
        rake_s = rake1
    else:
        stk_s = stk0
        dip_s = dip0
        rake_s = rake0
    return stk_s, dip_s, rake_s

# define a fault-normal vector
def faultnormalvec(stk, dip):
    nn = -np.sin(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    ne =  np.cos(np.deg2rad(stk)) * np.sin(np.deg2rad(dip))
    nd = -np.cos(np.deg2rad(dip))
    return np.array([ne, nn, nd])

from matplotlib.patches import Polygon
from matplotlib.animation import FuncAnimation


def fetch_JMA_aftershock():
    #URL = 'http://evrrss.eri.u-tokyo.ac.jp/cgi-bin/tsw-pro.pl?ACTION=LIST&DB=JMA1&PRMFILE=jma.deck.prm&TM1=2024/01/01&TM2=2024/01/07&RGN=125/25/148/48&DEP=0.0/700.0&MAG=0.0/9.9&SCALE=1.0&PROJECT=2&BORDER=1&GRID=0&CONF=&FMT=TXT'
    #response = requests.get(URL)
    #bs = BeautifulSoup(response.text, 'html.parser')
    #with open(os.path.join(datarootdir,'tseis.txt'), 'w') as file:
    #    pre=bs.find_all('pre')[0].text
    #    pre=pre[1:-1]
    #    file.write(pre)

    names = ['year','month','day','hour','minute','second','longitude','latitude','depth','magnitude']
    df_JMA = pd.read_csv(os.path.join(datarootdir,'tseis.txt'),header=None,names=names,delim_whitespace=True)
    df_JMA['datetime_JST'] = pd.to_datetime({
        'year': df_JMA['year'],
        'month': df_JMA['month'],
        'day': df_JMA['day'],
        'hour': df_JMA['hour'],
        'minute': df_JMA['minute'],
        'second': df_JMA['second'],
    })
    df_JMA['datetime_UTC'] = df_JMA['datetime_JST'] + pd.DateOffset(hours=-9)
    df_JMA
    df_JMA[ df_JMA['magnitude'] == np.max(df_JMA['magnitude']) ]
    return df_JMA

def fetch_JMA_background():
    #URL = 'http://evrrss.eri.u-tokyo.ac.jp/cgi-bin/tsw-pro.pl?ACTION=LIST&DB=JMA1&PRMFILE=jma.deck.prm&TM1=1900/01/01&TM2=2023/12/31&RGN=135/36/139/39&DEP=0.0/700.0&MAG=6.0/9.9&SCALE=1.0&PROJECT=2&BORDER=1&GRID=0&CONF=&FMT=TXT'
    #response = requests.get(URL)
    #bs = BeautifulSoup(response.text, 'html.parser')
    #with open(os.path.join(datarootdir,'tseis_background.txt'), 'w') as file:
    #    pre=bs.find_all('pre')[0].text
    #    pre=pre[1:-1]
    #    file.write(pre)

    names = ['year','month','day','hour','minute','second','longitude','latitude','depth','magnitude']
    df_JMA_bg = pd.read_csv(os.path.join(datarootdir,'tseis_background.txt'),header=None,names=names,delim_whitespace=True)
    df_JMA_bg['datetime_JST'] = pd.to_datetime({
        'year': df_JMA_bg['year'],
        'month': df_JMA_bg['month'],
        'day': df_JMA_bg['day'],
        'hour': df_JMA_bg['hour'],
        'minute': df_JMA_bg['minute'],
        'second': df_JMA_bg['second'],
    })
    df_JMA_bg['datetime_UTC'] = df_JMA_bg['datetime_JST'] + pd.DateOffset(hours=-9)
    return df_JMA_bg

def fetch_JMA_background_minor():
    #URL = 'http://evrrss.eri.u-tokyo.ac.jp/cgi-bin/tsw-pro.pl?ACTION=LIST&DB=JMA1&PRMFILE=jma.deck.prm&TM1=1900/01/01&TM2=2024/01/01&RGN=135/36/139/39&DEP=0.0/700.0&MAG=3.0/9.9&SCALE=1.0&PROJECT=2&BORDER=1&GRID=0&CONF=&FMT=TXT'
    #response = requests.get(URL)
    #bs = BeautifulSoup(response.text, 'html.parser')
    #with open(os.path.join(datarootdir,'tseis_background.txt'), 'w') as file:
    #    pre=bs.find_all('pre')[0].text
    #    pre=pre[1:-1]
    #    file.write(pre)

    names = ['year','month','day','hour','minute','second','longitude','latitude','depth','magnitude']
    df_JMA_bg = pd.read_csv(os.path.join(datarootdir,'tseis_background.txt'),header=None,names=names,delim_whitespace=True)
    df_JMA_bg['datetime_JST'] = pd.to_datetime({
        'year': df_JMA_bg['year'],
        'month': df_JMA_bg['month'],
        'day': df_JMA_bg['day'],
        'hour': df_JMA_bg['hour'],
        'minute': df_JMA_bg['minute'],
        'second': df_JMA_bg['second'],
    })
    df_JMA_bg['datetime_UTC'] = df_JMA_bg['datetime_JST'] + pd.DateOffset(hours=-9)
    return df_JMA_bg

def fetch_JMA_2007Noto():
    #URL = 'http://evrrss.eri.u-tokyo.ac.jp/cgi-bin/tsw-pro.pl?ACTION=LIST&DB=JMA1&PRMFILE=jma.deck.prm&TM1=2007/03/25&TM2=2007/03/26&RGN=135/36/139/39&DEP=0.0/700.0&MAG=0.0/9.9&SCALE=1.0&PROJECT=2&BORDER=1&GRID=0&CONF=&FMT=TXT'
    #response = requests.get(URL)
    #bs = BeautifulSoup(response.text, 'html.parser')
    #with open(os.path.join(datarootdir,'tseis_200Noto.txt'), 'w') as file:
    #    pre=bs.find_all('pre')[0].text
    #    pre=pre[1:-1]
    #    file.write(pre)

    names = ['year','month','day','hour','minute','second','longitude','latitude','depth','magnitude']
    df_JMA_bg = pd.read_csv(os.path.join(datarootdir,'tseis_200Noto.txt'),header=None,names=names,delim_whitespace=True)
    df_JMA_bg['datetime_JST'] = pd.to_datetime({
        'year': df_JMA_bg['year'],
        'month': df_JMA_bg['month'],
        'day': df_JMA_bg['day'],
        'hour': df_JMA_bg['hour'],
        'minute': df_JMA_bg['minute'],
        'second': df_JMA_bg['second'],
    })
    df_JMA_bg['datetime_UTC'] = df_JMA_bg['datetime_JST'] + pd.DateOffset(hours=-9)
    return df_JMA_bg

def fetch_JMA_2023Noto():
    #URL = 'http://evrrss.eri.u-tokyo.ac.jp/cgi-bin/tsw-pro.pl?ACTION=LIST&DB=JMA1&PRMFILE=jma.deck.prm&TM1=2020/11/01&TM2=2023/05/12&RGN=135/36/139/39&DEP=0.0/700.0&MAG=0.0/9.9&SCALE=1.0&PROJECT=2&BORDER=1&GRID=0&CONF=&FMT=TXT'
    #response = requests.get(URL)
    #bs = BeautifulSoup(response.text, 'html.parser')
    #with open(os.path.join(datarootdir,'tseis_2023Noto.txt'), 'w') as file:
    #    pre=bs.find_all('pre')[0].text
    #    pre=pre[1:-1]
    #    file.write(pre)

    names = ['year','month','day','hour','minute','second','longitude','latitude','depth','magnitude']
    df_JMA_bg = pd.read_csv(os.path.join(datarootdir,'tseis_2023Noto.txt'),header=None,names=names,delim_whitespace=True)
    df_JMA_bg['datetime_JST'] = pd.to_datetime({
        'year': df_JMA_bg['year'],
        'month': df_JMA_bg['month'],
        'day': df_JMA_bg['day'],
        'hour': df_JMA_bg['hour'],
        'minute': df_JMA_bg['minute'],
        'second': df_JMA_bg['second'],
    })
    df_JMA_bg['datetime_UTC'] = df_JMA_bg['datetime_JST'] + pd.DateOffset(hours=-9)
    return df_JMA_bg


def gen_shaped_text(n,m,tmpdf):
    text = str(n).rjust(4) + str(m).rjust(4) +\
    str('{:.10f}'.format(tmpdf['lat'].values[0])).rjust(17) + str('{:.10f}'.format(tmpdf['lon'].values[0])).rjust(17) +\
    str('{:.3f}'.format(tmpdf['dep'].values[0])).rjust(10) + str('{:.3f}'.format(tmpdf['strike'].values[0])).rjust(10) + str('{:.3f}'.format(tmpdf['dip'].values[0])).rjust(10) +\
    str('{:.3f}'.format(90)).rjust(10) + str('{:.2f}'.format(0)).rjust(10) + str('{:.3f}'.format(0)).rjust(10)
    return text


# scrape the coordinates of active fault segments
def read_kmz(fname):
    zfile = zipfile.ZipFile(fname)
    kml_string = zfile.read(zfile.filelist[0].filename)
    kml = KML()
    kml.from_string(kml_string)
    return kml

def load_active_faults():
    # https://unit.aist.go.jp/ievg/actfault-rg/kmz/newest/MG_trace_gbank_e_g.kmz
    kml = read_kmz('./MG_trace_gbank_e_g.kmz')
    coords_xys = []
    fault_ids = []
    for Document_feat in kml.features():  # 1 level.
        for Folder_feat in Document_feat.features():
            for SubFolder_feat in Folder_feat.features():  # 20 levels.
                for PlaceMark_feat in SubFolder_feat.features():
                    strings = PlaceMark_feat.geometry
                    # fault segment ID is something like; 183-06 and can be found at e.g.)
                    # https://gbank.gsj.jp/activefault/segment_param_e?SearchTYPE=&fval_type1=183-06&segment_id=183-06&topic_list=2&search_mode=2
                    # you can check the detail of fault segment by; print(PlaceMark_feat.description)
                    #if PlaceMark_feat.name == '183-06':
                    fault_ids.append(PlaceMark_feat.name)
                    for string in strings.geoms:
                        #coords_xys.append(np.vstack(string.coords.xy))
                        coords_xys.append(np.vstack(string.coords))
    #print('Total number of fault segments:', len(coords_xys))
    return coords_xys

from cmaptools import readcpt, joincmap, DynamicColormap


def get_plunge_azimuth(focmec):
    '''
    input: 6 components of  mooment tensor in USE convention (e.g., default GCMT, psmeca format)
    output: A list of [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth] angles in degree
    '''
    m1,m2,m3,m4,m5,m6 = focmec[0],focmec[1],focmec[2],focmec[3],focmec[4],focmec[5]
    mt = obspy.imaging.beachball.MomentTensor(m1,m2,m3,m4,m5,m6, 26)

    (d, v) = np.linalg.eigh(mt.mt)
    pl = np.arcsin(-v[0])
    az = np.arctan2(v[2], -v[1])
    for i in range(0, 3):
        if pl[i] <= 0:
            pl[i] = -pl[i]
            az[i] += np.pi
        if az[i] < 0:
            az[i] += 2 * np.pi
        if az[i] > 2 * np.pi:
            az[i] -= 2 * np.pi

    pl = np.rad2deg(pl)
    az = np.rad2deg(az)

    p_plunge = pl[0]
    p_azimuth = az[0]
    n_plunge = pl[1]
    n_azimuth = az[1]
    t_plunge = pl[2]
    t_azimuth = az[2]

    return [p_plunge, p_azimuth, n_plunge, n_azimuth, t_plunge, t_azimuth]

def load_coastline(resolution='i'):
    # availabla via https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/
    xs_c,ys_c = [],[]
    for num in np.arange(1,7,1):
        src = shapefile.Reader('/Users/onion/Downloads/gshhg-shp-2.3.7/GSHHS_shp/'+resolution+'/GSHHS_'+resolution+'_L'+str(num)+'.shp')
        for tmp in src.shapeRecords():
            x, y = [i[0] for i in tmp.shape.points[:]], [i[1] for i in tmp.shape.points[:]]
            if np.min(x) >= 130 and np.max(x) <= 150 and np.min(y) >= 30 and np.max(y) <= 42:
                xs_c.append(x)
                ys_c.append(y)
    return xs_c,ys_c

