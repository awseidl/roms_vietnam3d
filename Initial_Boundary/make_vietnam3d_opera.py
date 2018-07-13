import netCDF4
from datetime import datetime, timedelta
import ROMStools
import numpy as np
from numpy import f2py
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import depths
import time
from calendar import monthrange
import argparse
import sys
#import matplotlib.pyplot as plt

def get_var(nc,var,index):
    nc = netCDF4.Dataset(nc)
    for i in range(len(var)):
        if var[i] in nc.variables.keys():
            if ( len(nc.variables[var[i]].shape) == 3 ):
                yd, xd = nc.variables[var[i]][0,:,:].shape
                varo = np.zeros([1,yd,xd])
                varo[0,:,:]=nc.variables[var[i]][index,:,:]
            elif ( len(nc.variables[var[i]].shape) == 4 ):
                sd, yd, xd = nc.variables[var[i]][0,:,:,:].shape
                varo = np.zeros([1,sd,yd,xd])
                varo[0,:,:,:] = nc.variables[var[i]][index,:,:,:]
    nc.close()
    return varo

def fill(var):
    fillvalue = var.min()
    varmasked = np.ma.masked_values(var,fillvalue)
    if ( len(var.shape) == 4 ):
        for s in range(len(var[0,:,0,0])):
            var[0,s,:,:] = np.where((var[0,s,:,:] == fillvalue), varmasked[0,s,:,:].mean(), var[0,s,:,:])
            var[0,s,:,:] = gaussian_filter(var[0,s,:,:],sigma=1)
    elif ( len(var.shape) == 3 ):
        var[0,:,:] = np.where((var[0,:,:] == fillvalue), varmasked[0,:,:].mean(), var[0,:,:])
        var[0,:,:] = gaussian_filter(var[0,:,:],sigma=1)
    return var

def fill_xdir(var):
    fillvalue = var.min()
    varmasked = np.ma.masked_values(var,fillvalue)
    if ( len(var.shape) == 4 ):
        for s in range(len(var[0,:,0,0])):
            for j in range(len(var[0,s,:,0])):
                var[0,s,j,:] = np.where((var[0,s,j,:] == fillvalue), varmasked[0,s,j,:].mean(), var[0,s,j,:])
            var[0,s,:,:] = gaussian_filter(var[0,s,:,:],sigma=1)
    elif ( len(var.shape) == 3 ):
        for j in range(len(var[0,:,0])):
            var[0,j,:] = np.where((var[0,j,:] == fillvalue), varmasked[0,j,:].mean(), var[0,j,:])
        var[0,:,:] = gaussian_filter(var[0,:,:],sigma=2)
    return var

### MAIN ####
#yyyy,mm,hh = read_input()

daysoffset = 0
ifile = sys.argv[1]
gfile = 'east_sea_grd-smooth.nc'
#afile = 'http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m-anglematrix/angle_norkyst-800m_grd.nc'
odir  = 'Output'

#Output grid parameters:
N = 20
hc = 100
theta_b = 0.3
theta_s = 6.0
Vtransform = 2
Vstretching = 4

latro, lonro = ROMStools.read_grid_info(gfile)
latri, lonri = ROMStools.read_grid_info(ifile)

#fig, axs = plt.subplots(1,2)
#grdin = axs[0].contourf(latri)
#fig.colorbar(grdin, ax=axs[0])
#grdout = axs[1].contourf(latro)
#fig.colorbar(grdout, ax=axs[1])
#plt.show()
#exit()

#anglei = ROMStools.get_var_pointer(afile, ['angle'])[1:-1,1:-1]
#anglei = ROMStools.get_var_pointer(afile, ['angle'])[:,:]
angleo = ROMStools.get_var_pointer(gfile, ['angle'])[:,:]

h = ROMStools.get_var_pointer(gfile, 'h')

timee = ROMStools.get_var_pointer(ifile,['ocean_time','time'])
time2 = netCDF4.num2date(timee[:],timee.units)
idate = datetime.now()-timedelta(days=daysoffset)
#starttime = ROMStools.find_nearest_index(timee[:], netCDF4.date2num(idate,units=timee.units,calendar=timee.calendar))
starttime = 0
for ii in range(starttime,len(timee[:])):
    try:
        ofile    = odir+'/vietnam3d_%02d.nc' % (ii+1)
	print ofile
        print 'Date: '+str(time2[ii])

        # Get variables
        zetai = get_var(ifile,['zeta','SSH'],ii)
        salti = get_var(ifile,['salinity','salt'],ii)
        tempi = get_var(ifile,['temperature','temp'],ii)
        ui = get_var(ifile,['u'],ii)
        vi = get_var(ifile,['v'],ii)
        ubari = get_var(ifile,['ubar'],ii)
        vbari = get_var(ifile,['vbar'],ii)

        # Fill masked values
        ubari = fill(ubari)
        vbari = fill(vbari)
        ui = fill(ui)
        vi = fill(vi)
        #zetai = fill_xdir(zetai)
        salti = fill(salti)
        tempi = fill(tempi)

        #print '----------------------'
        #print 'vector rotation'
        uie = ui
	vin = vi
        ubarie = ubari
	vbarin = vbari
        #print 'done...'
        #print '----------------------'


        print '----------------------'
        print 'horizontal interpolation'
        zetao_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, zetai)
        zetao_tmp = fill_xdir(zetao_tmp)
        salto_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, salti)
        tempo_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, tempi)
        uo_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, uie) # u on thredds is not staggered
        vo_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, vin) # v on thredds is not staggered
        ubaro_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, ubarie) # ubar on thredds is not staggered
        #ubaro_tmp = fill_xdir(ubaro_tmp)
        vbaro_tmp = ROMStools.hor_interp(latri, lonri, latro, lonro, vbarin) # vbar on thredds is not staggered
        #vbaro_tmp = fill_xdir(vbaro_tmp)
        print 'done...'
        print '----------------------'

        #Delete unused arrays:
        del(zetai)
        del(salti)
        del(ui)
        del(ubari)
        del(vi)
        del(vbari)

        zetao = zetao_tmp
        ubaro = ubaro_tmp
        vbaro = vbaro_tmp

        print '----------------------'
        print 'vector rotation'
        uo_tmp_g, vo_tmp_g = ROMStools.rotate_vectors_togrid(angleo,uo_tmp,vo_tmp)
        ubaro_g, vbaro_g = ROMStools.rotate_vectors_togrid(angleo,ubaro,vbaro)
        print 'done...'
        print '----------------------'
        #Delete unused arrays:
        del(ubaro_tmp)
        del(vbaro_tmp)

        print '----------------------'
        print 'vertical interpolation'
        print '----------------------'
        salto=np.zeros([len(zetao[:,0,0]),N,len(h[:,0]),len(h[0,:])])
        tempo=np.zeros([len(zetao[:,0,0]),N,len(h[:,0]),len(h[0,:])])
        uo=np.zeros([len(zetao[:,0,0]),N,len(h[:,0]),len(h[0,:])])
        vo=np.zeros([len(zetao[:,0,0]),N,len(h[:,0]),len(h[0,:])])

        for t in range(len(zetao[:,0,0])): #Loop over time
            depth3di = ROMStools.read_vert_grid_from_file(ifile,h)
            depth3do = -depths.get_zrho(Vtransform, Vstretching, N, theta_s, theta_b, h[:,:], hc, zetao_tmp[t,:].squeeze())
            salto[t,:] = ROMStools.vert_int3d(depth3di, depth3do, salto_tmp[t,:].squeeze())
            tempo[t,:] = ROMStools.vert_int3d(depth3di, depth3do, tempo_tmp[t,:].squeeze())
            uo[t,:] = ROMStools.vert_int3d(depth3di, depth3do, uo_tmp_g[t,:,:,:].squeeze()) #hack!
            vo[t,:] = ROMStools.vert_int3d(depth3di, depth3do, vo_tmp_g[t,:,:,:].squeeze()) #hack!
        print 'done...'
        print '----------------------'

        print '-------------------'
        print 'Writing curvilinear file...'
        ROMStools.create_climfile(ofile,['salt','temp','zeta','u','v','ubar','vbar'],
                                     [salto,tempo,zetao,0.5*(uo[:,:,:,:-1]+uo[:,:,:,1:]),0.5*(vo[:,:,:-1,:]+vo[:,:,1:,:]),
                                      0.5*(ubaro_g[:,:,:-1]+ubaro_g[:,:,1:]),0.5*(vbaro_g[:,:-1,:]+vbaro_g[:,1:,:])])
        # Add clim_time variable to climfile:
        print 'curvilinear file ok'
        r = netCDF4.Dataset(ofile,'r+')
        vvar = r.createVariable('ocean_time','double',('ocean_time'))
        vvar[:] = timee[ii]
        r.close()
        print 'soon...'
        # def add_attribute_to_fil
        ROMStools.get_var_pointer(ofile,['ocean_time'],readwrite='r+').units = timee.units
        print '----------------------'
        #Delete more:
        del(tempo)
        del(salto)
        del(zetao)
        del(uo)
        del(vo)
        del(ubaro_g)
        del(vbaro_g)
        print 'done'
    except Exception as ex:
        print 'Error: '+str(ex)


