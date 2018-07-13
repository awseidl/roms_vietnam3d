import numpy as np
import math
from datetime import datetime,timedelta
import multiprocessing
import netCDF4
import sys
#import scipy
from scipy.ndimage import gaussian_filter
from scipy.interpolate import griddata
import os
'''
Module containing tools for ROMS.
nilsmk@met.no - 2014
'''

def vert_int3d(z3di, z3do, var3di):
    '''
    Function for interpolating 3D field between depths given
    in two 3D matrixes.
    '''
    yo, xo = var3di[0,:,:].shape

    # Remove nans:
    idx = np.isfinite(var3di)
    varmin = np.nanmin(var3di)
    varmax = np.nanmax(var3di)
    if var3di[1,0,0] < var3di[0,0,0]:
        var3di = np.where(idx,var3di,varmin)
    else:
        var3di = np.where(idx,var3di,varmax)
    #z3di = np.where(idx,z3di,z3di)

    # Test of depth values are pos or neg:
    if z3di[0,0,1] < z3di[0,0,0]:
        print 'Swapping depths...'
        var3di = var3di[::-1,...]

    if len(z3do.shape) == 3:
        var3do = np.zeros_like(z3do)
    elif len(z3do.shape) == 1:
        var3do = np.zeros([len(z3do), yo, xo])
    for j in range(yo):
        for i in range(xo):
            if len(z3do.shape) == 3:
                var3do[:,j,i] = np.interp(z3do[:,j,i], z3di[:,j,i], var3di[:,j,i])
            elif len(z3do.shape) == 1:
                var3do[:,j,i] = np.interp(z3do[:], z3di[:,j,i], var3di[:,j,i])
    return var3do

def getKE(nc,grd,timevar,nlev):
    '''
    Kinetic energy timeseries from 3D array.
    '''
    mask = grd.variables['mask_rho'][:,:]
    ocean_time = nc.variables[timevar][:]
    KESt = np.zeros(len(ocean_time))
    print "Looping over time..."
    for t in range(len(ocean_time)):
        if (nlev > 100):
            print "error... program not ready for this yet..."
            exit
        else:
            u_staggS = (nc.variables['u'][t,nlev,:,:])*100 # change to cms-1
            uvarS = ((u_staggS[0:-1,:] + u_staggS[1:,:]) / 2) * mask[1:,1:]
            v_staggS = (nc.variables['v'][t,nlev,:,:])*100 # change to cms-1
            vvarS = ((v_staggS[:,0:-1] + v_staggS[:,1:]) / 2) * mask[1:,1:]
        KES = 0.5 * ((uvarS**2) + (vvarS**2))
        KESt[t] = np.mean(KES[:,:])
    return KESt


def getEKE(nc,grd,timevar,nlev,avgperiod):
    '''
    Eddy kinetic energy timeseries from 3D array.
    '''
    mask = grd.variables['mask_rho'][:,:]
    ocean_time = nc.variables[timevar][:]
    EKESt = np.zeros(len(ocean_time))
    MKESt = np.zeros(len(ocean_time))
    print "Looping over time... Running mean of "+str(avgperiod)+" days"
    for t in range(len(ocean_time)):
        if (nlev > 100):
            print "error... program not ready for this yet..."
            exit
        else:
            if (t < avgperiod/2):
                #print "opt1, "+str(t)
                u_staggS = (nc.variables['u'][0:(t+(avgperiod/2)),nlev,:,:])*100 # change to cms-1
                v_staggS = (nc.variables['v'][0:(t+(avgperiod/2)),nlev,:,:])*100 # change to cms-1
            elif ((len(ocean_time)-t) < avgperiod/2):
                #print "opt2, "+str(t)
                u_staggS = (nc.variables['u'][(t-(avgperiod/2)):,nlev,:,:])*100 # change to cms-1
                v_staggS = (nc.variables['v'][(t-(avgperiod/2)):,nlev,:,:])*100 # change to cms-1
            else:
                #print "opt3, "+str(t)
                u_staggS = (nc.variables['u'][(t-(avgperiod/2)):(t+(avgperiod/2)),nlev,:,:])*100 # change to cms-1
                v_staggS = (nc.variables['v'][(t-(avgperiod/2)):(t+(avgperiod/2)),nlev,:,:])*100 # change to cms-1
            uvarS = ((u_staggS[:,0:-1,:] + u_staggS[:,1:,:]) / 2)
            vvarS = ((v_staggS[:,:,0:-1] + v_staggS[:,:,1:]) / 2)
            umeanS = np.mean(uvarS,axis=0) * mask[1:,1:]
            vmeanS = np.mean(vvarS,axis=0) * mask[1:,1:]
        KES = 0.5 * np.mean((uvarS**2) + (vvarS**2),axis=0) * mask[1:,1:]
        MKES = 0.5 * ((umeanS**2) + (vmeanS**2))
        EKESt[t] = np.mean(KES[:,:] - MKES)
        MKESt[t] = np.mean(MKES)
    return EKESt,MKESt

def calc_energy_2d(u, v):
    umean = np.mean(u,axis=0)
    vmean = np.mean(v,axis=0)
    KE = 0.5 * np.mean((u**2 + v**2),axis=0).squeeze()
    MKE = 0.5 * ((umean**2) + (vmean**2))
    EKE = KE - MKE
    return KE, MKE, EKE

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def destagger_uv(u, v, scale=1, nlev=-1):
    '''
    Function to destagger u and v
    If nlev == -2, assume input field is 2D (e.g. ubar/vbar)
    If nlev >=  0, return a vertical layer from 3D field as a 2D field
    If nlev == -3, return 3D matrix
    '''
    if len(u.shape) == 3:
        ud = ((u[:,1:-1,0:-1] + u[:,1:-1,1:]) / 2)*scale
        vd = ((v[:,0:-1,1:-1] + v[:,1:,1:-1]) / 2)*scale
    elif len(u.shape) == 4:
        ud = ((u[:,nlev,1:-1,0:-1] + u[:,nlev,1:-1,1:]) / 2)*scale
        vd = ((v[:,nlev,0:-1,1:-1] + v[:,nlev,1:,1:-1]) / 2)*scale
    return ud,vd

def destagger_uv2d(u, v):
    '''
    Function to destagger u and v
    '''
    ud = ((u[1:-1,0:-1] + u[1:-1,1:]) / 2)
    vd = ((v[0:-1,1:-1] + v[1:,1:-1]) / 2)
    return ud,vd

def rmse(A,B):
    rmse = 0.0
    for i in range(len(A)):
        rmse = rmse + ((A[i]-B[i])**2)
    rmse = math.sqrt(rmse/len(A))
    return rmse

def get_tseries_2Dvar(nc,variable,timevar,x,y):
    Tseries = nc.variables[variable][:,y,x]
    timevar2D = nc.variables[timevar][:]
    return Tseries, timevar2D

def get_tseries_3Dvar(nc,variable,timevar,x,y,s=-1):
    Tseries = nc.variables[variable][:,s,y,x]
    timevar2D = nc.variables[timevar][:]
    return Tseries, timevar2D

def sec2str(secarray):
    D1970 = datetime.strptime('1970-01-01','%Y-%m-%d')
    outtime=[(D1970 + timedelta(seconds=s)) for s in secarray]
    return outtime

def calc_mean_big_array3D(nc,variable,nlev,start,stopp):
    n = 0
    meanarray = nc.variables[variable][start,nlev,:,:].squeeze()
    for i in range(start+1,stopp+1):
        print i
        meanarray = meanarray + nc.variables[variable][i,nlev,:,:].squeeze()
        n = n +1
    meanarray = meanarray / n
    return meanarray

def calc_mean_big_array2D(nc,variable,start,stopp):
    n = 0
    meanarray = nc.variables[variable][start,:,:].squeeze()
    for i in range(start+1,stopp+1):
        #print i
        meanarray = meanarray + nc.variables[variable][i,:,:].squeeze()
        n = n +1
    meanarray = meanarray / n
    return meanarray

# def calc_meanKE_big_array3D(nc,nlev,start,stopp):
#     u, v = destagger_uv(nc.variables['u'], nc.variables['v'], scale=1, nlev)
#     n = 0
#     meanarray = ((u[start,nlev,:,:].squeeze()**2 + v[start,nlev,:,:].squeeze()**2) / 2)
#     for i in range(start+1,stopp+1):
#         #print i
#         meanarray = meanarray + (u[i,nlev,:,:].squeeze())
#         n = n +1
#     meanarray = meanarray / n
#     return meanarray
#

###################################################################################################
# The following stuff was developed for generating clim from thredds:

#Define some functions, put these into somewhere else?? ROMStools maybe?
# def get_var(nc,var):
#     #General get_var-function that will return variable of any dimension
#     nc = netCDF4.Dataset(nc)
#     for i in range(len(var)):
#         if var[i] in nc.variables.keys():
#             print 'read '+var[i]
#             varo = nc.variables[var[i]][:]
#             break
#     return varo


def get_var_pointer(fname,var,readwrite='r'):
    #Returns a pointer to a variable
    nc = netCDF4.Dataset(fname,readwrite)
    foundvar = False
    for i in range(len(var)):
        if var[i] in nc.variables.keys():
            print 'read '+var[i]+' as pointer...'
            varo = nc.variables[var[i]]
            foundvar = True
            break
    if ( foundvar == False ):
        print 'Var '+var+' not found on file, exit'
        sys.exit()
    #nc.close()
    return varo

def get_var(nc,var,start=0,stopp=0,skiptime=0):
    #Function used for 4D-variables, i.e. [time,depth,y,x]
    if ((stopp - start) < 0):
        print 'Stopping time less than starting time, exit'
        sys.exit()
    else:
        nc = netCDF4.Dataset(nc)
        for i in range(len(var)):
            if var[i] in nc.variables.keys():
                print 'read '+var[i]+', dims = '+str(len(nc.variables[var[i]].shape))
                if ( len(nc.variables[var[i]].shape) == 4 ):
                    if ((stopp - start) >= 1):
                        # Read variable as time chuncks:
                        sd, yd, xd = nc.variables[var[i]][0,:,:,:].shape
                        varo = np.zeros([stopp-start+1,sd,yd,xd])
                        for t in range((stopp-start+1)):
                            varo[t,:,:,:] = nc.variables[var[i]][start+t,:,:,:]
                    else:
                        print 'problem'
                        sys.exit()
                        #varo = nc.variables[var[i]][start:stopp,:,:,:]
                elif ( len(nc.variables[var[i]].shape) == 3 ):
                    if ((stopp - start) >= 1):
                        # Read variable as time chuncks:
                        yd, xd = nc.variables[var[i]][0,:,:].shape
                        varo = np.zeros([stopp-start+1,yd,xd])
                        for t in range((stopp-start+1)):
                            varo[t,:,:] = nc.variables[var[i]][start+t,:,:]
                    else:
                        print 'problem'
                        sys.exit()
                        #varo = nc.variables[var[i]][start:stopp,:,:]
                elif ( len(nc.variables[var[i]].shape) == 2 ):
                    varo = nc.variables[var[i]][:,:]
                elif ( len(nc.variables[var[i]].shape) == 1 ):
                    varo = nc.variables[var[i]][:]
                    print 'expanding with meshgrid...'
                else:
                    print 'error on dimensions, exit'
                    sys.exit()
                break
    return varo

def read_grid_info(grdname):
    #will read all info from a gridfile/romsfile
    latr = get_var(grdname,['lat_rho','latitude','ulat'])
    lonr = get_var(grdname,['lon_rho','longitude','ulon'])
    if ( len(latr.shape) == 1 or len(lonr.shape) == 1 ):
        lonr2,latr2 = np.meshgrid(lonr,latr)
        return latr2, lonr2
    else:
        return latr, lonr

def make_bry(clmfile):
    #will make bryfile from climfile
    return

def hor_interp(lati,loni,lato,lono,vari,method='nearest'):
    #This will interpolate the input variable to the output grid
    y, x = lato.shape
    #method='linear'
    if ( len(vari.shape) == 2 ):
        #vari = fill(vari)
        varo = griddata((np.hstack(lati),np.hstack(loni)),np.hstack(vari),(lato,lono), method)
    elif ( len(vari.shape) == 3 ):
        t, ydum, xdum = vari.shape
        varo = np.zeros([t,y,x])
        for i in range(t):
            #vari[i,:,:] = fill(vari[i,:,:])
            varo[i,:,:] = griddata((np.hstack(lati),np.hstack(loni)),np.hstack(vari[i,:,:]),(lato,lono), method)
    elif ( len(vari.shape) == 4 ):
        t, s, ydum, xdum = vari.shape
        varo = np.zeros([t,s,y,x])
        for i in range(t):
            for j in range(s):
                #vari[i,j,:,:] = fill(vari[i,j,:,:])
                varo[i,j,:,:] = griddata((np.hstack(lati),np.hstack(loni)),np.hstack(vari[i,j,:,:]),(lato,lono), method)
    else:
        print 'unsupported number of dims on variable '+str(len(vari.shape))
        sys.exit()
    return varo

# def vert_interp(var,vert_coord=Z_LEVELS):
#     if vert_coord == Z_LEVELS:
#     #Vertical interpolation from z-levels to S-levels:
#     #Should read input depths from ifile..
#         dosomething = 0
#     elif vert_coord == S_LEVELS:
#     #Vertical interpolation from S-levels to S-levels:
#     #Where should this info be read from? could do similar to roms_ini_from_profile, but then needs octant...
#         dosomething = 1
#     else:
#         print 'Illegal vertical coordinate, use either Z_DEPTH or S_LEVELS'
#         sys.exit()
#     return var

def filterfunc(var2d, criteria=0.5, glattvar=8):
    # Maybe crit and glattvar should be function of max/man gradients?
    while (np.max(np.gradient(var2d)) > criteria):
        #print 'filtering...'
        var2d = gaussian_filter(var2d, sigma=glattvar)
    return var2d

def fill(var, fillvalue, criteria=0.5, glattvar=8):
    import inpaint
    varmasked = np.ma.masked_less(var,fillvalue)
    if ( len(var.shape) == 4 ):
        #Depths should be 2nd dimension...
        #print 'looping over time and depths'
        for t in range(len(var[:,0,0,0])):
            #for s in reversed(range(len(var[0,:,0,0])-1)):
            for s in range(len(var[0,:,0,0])):
                var[t,s,:,:] = np.where((var[t,s,:,:] <= fillvalue), varmasked[t,s,:,:].mean(), var[t,s,:,:])
                #print 'Gradients before: '+str(np.max(np.gradient(var[t,s,:,:])))+', '+str(np.min(np.gradient(var[t,s,:,:])))
                var[t,s,:,:] = filterfunc(var[t,s,:,:], criteria=0.05)
                #print 'Gradients after: '+str(np.max(np.gradient(var[t,s,:,:])))+', '+str(np.min(np.gradient(var[t,s,:,:])))
                if (s<(len(var[0,:,0,0])-1)): # this should fill values below seabed, but doesn't seem to work 100%
                    #print s
                    var[t,s+1,:,:] = np.where((var[t,s+1,:,:] <= fillvalue), var[t,s,:,:], var[t,s+1,:,:])
    elif ( len(var.shape) == 3 ):
        #print 'looping over time'
        for t in range(len(var[:,0,0])):
            var[t,:,:] = np.where(var[t,:,:] <= fillvalue, varmasked.mean(), var[t,:,:])
            #print 'Gradients before: '+str(np.max(np.gradient(var[t,:,:])))+', '+str(np.min(np.gradient(var[t,:,:])))
            var[t,:,:] = filterfunc(var[t,:,:], criteria=0.01)
            #print 'Gradients after: '+str(np.max(np.gradient(var[t,:,:])))+', '+str(np.min(np.gradient(var[t,:,:])))
    else:
        print 'error!!!'
    return var

def create_climfile(ofile,varnames,varlist):
    if ( len(varnames) != len(varlist) ):
        print 'check list of varnames and vars'
        sys.exit()
    # Save data to netCDF-file:
    try:
        if os.path.exists(ofile):
            print 'Outfile exists, will delete it'
            os.remove(ofile)
    #             print 'Outfile exists, will exit'
    #             sys.exit()
        r = netCDF4.Dataset(ofile,'w',type='NETCDF4-CLASSIC')
        for i in range(len(varlist)):
            if ( i == 0 ):
                #Add dimensions
                tdim, sdim, ydim, xdim = varlist[i].shape

                # 4 dimensions, time is unlimited
                dt  = r.createDimension('ocean_time',None)
                dxr = r.createDimension('xi_rho',xdim)
                dyr = r.createDimension('eta_rho',ydim)
                dxr = r.createDimension('xi_u',xdim-1)
                dyr = r.createDimension('eta_v',ydim-1)
                dsr = r.createDimension('s_rho',sdim)

                #vdt = r.createVariable('time','double',('time'))
                #vdt.units='seconds since 1970-01-01 00:00:00' #time.units
                # one global attribute
                r.title = "Made by andrewws@met.no with assistance from nilsmk@met.no"
            #Add variable
            add_var_to_ncobj(r,varnames[i],varlist[i])
            r.sync()
        r.close()
    #             else:
    #                 print 'error with dimensions, exit'
    #                 sys.exit()
    except Exception as ex:
        print "Could not prepare netCDF file for writing (%s)" % (ex,)
    return

def add_var_to_ncobj(r,varname,var):
    print len(var.shape)
    print 'varname: '+varname
    if ( len(var.shape) == 1 ):
        print 'Assume time dim'
        if ( varname == 'ocean_time' ):
            vvar = r.createVariable('ocean_time','double',('ocean_time'))
        elif ( varname == 's_rho' ):
            vvar = r.createVariable('s_rho','double',('s_rho'))
    elif ( len(var.shape) == 2 ):
        print 'Assume x,y dims'
        if ( varname == 'u' or varname == 'ubar'):
            vvar = r.createVariable(varname,'d',('eta_rho','xi_u'),fill_value=1e+37)
        elif ( varname == 'v' or varname == 'vbar'):
            vvar = r.createVariable(varname,'d',('eta_v','xi_rho'),fill_value=1e+37)
        else:
            vvar = r.createVariable(varname,'d',('eta_rho','xi_rho'),fill_value=1e+37)
    elif ( len(var.shape) == 3 ):
        print 'Assume x,y,time dims'
        if ( varname == 'u' or varname == 'ubar'):
            vvar = r.createVariable(varname,'d',('ocean_time','eta_rho','xi_u'),fill_value=1e+37)
        elif ( varname == 'v' or varname == 'vbar'):
            vvar = r.createVariable(varname,'d',('ocean_time','eta_v','xi_rho'),fill_value=1e+37)
        else:
            vvar = r.createVariable(varname,'d',('ocean_time','eta_rho','xi_rho'),fill_value=1e+37)
    elif ( len(var.shape) == 4 ):
        print 'Assume x,y,s,time dims'
        if ( varname == 'u' or varname == 'ubar'):
            vvar = r.createVariable(varname,'d',('ocean_time','s_rho','eta_rho','xi_u'),fill_value=1e+37)
        elif ( varname == 'v' or varname == 'vbar'):
            vvar = r.createVariable(varname,'d',('ocean_time','s_rho','eta_v','xi_rho'),fill_value=1e+37)
        else:
            vvar = r.createVariable(varname,'d',('ocean_time','s_rho','eta_rho','xi_rho'),fill_value=1e+37)
    else:
        print 'error!'
    vvar[:] = var

def add_var_to_file(file,varname,var):
    r = netCDF4.Dataset(file,'r+')
    if ( len(var.shape) == 1 ):
        print 'Assume time dim'
        if ( varname == 'ocean_time' ):
            vvar = r.createVariable('ocean_time','double',('ocean_time'))
    vvar[:] = var

# def add_attribute_to_file(file,varname,attname,att):
#     r = netCDF4.Dataset(file)
#     if ( varname == 'ocean_time' ):
#         var = r.variables[varname]


def read_vert_grid_from_file(filename,h=None):
    # Will test if vertical levels are in meters, or s-levels??
    depth = get_var_pointer(filename, ['depth','s_rho'] )
    #print depth
    if h == None:
        y, x = get_var_pointer(filename, 'h')[:,:].shape
    else:
        y, x = h.shape
    print x,y
    # Fill the arrays:
    if ( depth.standard_name == 'depth' ):
        depth3d = np.zeros([len(depth),y,x])
        for i in range(len(depth)):
            depth3d[i,:,:] = depth[i]
    elif ( depth.standard_name == 'ocean_s_coordinate_g2' ):
        #print 'slevels not ready yet, exit'
        #sys.exit()
        #Get stuff needed for vertical interpolation:
        nc = netCDF4.Dataset(filename)
        cs_r = nc.variables['Cs_r'][:]
        hc = nc.variables['hc'][:]
        sc_r = nc.variables['s_rho'][:]
        sdim = len(sc_r)
        Vtrans = nc.variables['Vtransform'][:]
        depth3d = np.zeros((sdim,y,x))
        if (Vtrans == 1):
            for k in range(0,sdim):
                depth3d[k,:,:]   = hc * sc_r[k] + (h-hc) * cs_r[k]
        elif (Vtrans == 2):
            for k in range(0,sdim):
                depth3d[k,:,:]   = h * ((hc * sc_r[k] + h * cs_r[k]) / (h+hc) )
        else:
            print 'Unknown Vtrans - exit'
            sys.exit()
    else:
        print 'unknown vertical coord, exit'
        sys.exit()
    return depth3d

def calc_vert_grid(Vtransform, Vstretching, N, theta_s, theta_b, h, hc, zeta):
    pass

def rotate_vectors_tonorth(angle,ui,vi):
    #Rotate vectors
    cosa = np.cos(angle)
    sina = np.sin(angle)
    uo = np.zeros_like(ui)
    vo = np.zeros_like(vi)
    if ( len(ui.shape) == 4 and len(vi.shape) == 4 and ui.shape == vi.shape ):
        #Depths should be 2nd dimension...
        print 'looping over time and depths'
        for t in range(len(ui[:,0,0,0])):
            for s in range(len(ui[0,:,0,0])):
                uo[t,s,:,:] = (ui[t,s,:,:] * cosa) - (vi[t,s,:,:] * sina)
                vo[t,s,:,:] = (ui[t,s,:,:] * sina) + (vi[t,s,:,:] * cosa)
    elif ( len(ui.shape) == 3 and len(vi.shape) == 3 and ui.shape == vi.shape ):
        #Depths should be 2nd dimension...
        print 'looping over time'
        for t in range(len(ui[:,0,0])):
            uo[t,:,:] = (ui[t,:,:] * cosa) - (vi[t,:,:] * sina)
            vo[t,:,:] = (ui[t,:,:] * sina) + (vi[t,:,:] * cosa)
    elif ( len(ui.shape) == 2 and len(vi.shape) == 2 and ui.shape == vi.shape ):
        #Depths should be 2nd dimension...
        print 'no loops'
        uo[:,:] = (ui[:,:] * cosa) - (vi[:,:] * sina)
        vo[:,:] = (ui[:,:] * sina) + (vi[:,:] * cosa)
    else:
        print 'error!!!'
    return uo,vo

def rotate_vectors_togrid(angle,ui,vi):
    #Rotate vectors
    uo, vo = rotate_vectors_tonorth(-angle,ui,vi)
    return uo,vo

def flat_fill(varin, crit=-10000000000):
    if crit == -10000000000:
        crit = varin.min()
    mask=(varin==crit)
    varin[mask]=np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), varin[~mask])
    return varin

def movingaverage(interval, window_size):
    #Stolen from http://stackoverflow.com/questions/11352047/finding-moving-average-from-data-points-in-python
    window = np.ones(int(window_size))/float(window_size)
    glidemiddel=np.convolve(interval, window, 'same')
    glidemiddel[:int(window_size/2)] = interval[:int(window_size/2)]
    glidemiddel[-int(window_size/2):] = interval[-int(window_size/2):]
    return glidemiddel

def calc_volume_ROMS(nc,y1=1,y2=-1,x1=1,x2=-1):
    #Get stuff needed for vertical interpolation:
    H = nc.variables['h'][y1:y2,x1:x2]
    pm = nc.variables['pm'][y1:y2,x1:x2]
    pn = nc.variables['pn'][y1:y2,x1:x2]
    dx = 1 / pm
    dy = 1 / pn
    V = np.sum(dx * dy * H)
    return V

def calc_area_ROMS(nc,y1=1,y2=-1,x1=1,x2=-1):
    #Get stuff needed for vertical interpolation:
    pm = nc.variables['pm'][y1:y2,x1:x2]
    pn = nc.variables['pn'][y1:y2,x1:x2]
    dx = 1 / pm
    dy = 1 / pn
    A = np.sum(dx * dy)
    return A

def return_grid_dims_ROMS(nc,y1=1,y2=-1,x1=1,x2=-1):
    #Get stuff needed for vertical interpolation:
    pm = nc.variables['pm'][y1:y2,x1:x2]
    pn = nc.variables['pn'][y1:y2,x1:x2]
    dx = 1 / pm
    dy = 1 / pn
    return dx,dy

def qdrag_calc_ROMS(nc,tval,s=0,rdrg2=0.003):
    '''
    Code is partly/mostly stolen from ROMS/Nonlinear/set_vbc.F
    '''
    u = nc.variables['u'][tval,s,:,:]
    v = nc.variables['v'][tval,s,:,:]
    yudim, xudim = u.shape
    yvdim, xvdim = v.shape
    bustr = np.zeros([yudim,xudim])
    bvstr = np.zeros([yvdim,xvdim])
    #
    cff1 = 0.25 * (v[:yudim-2,1:xudim]+v[1:yudim-1,1:xudim]+v[:yudim-2,:xudim-1]+v[1:yudim-1,:xudim-1])
    cff2 = np.sqrt((u[:yudim-2,1:xudim]**2)+(cff1**2))
    bustr[:yudim-2,1:xudim] = rdrg2 * u[:yudim-2,1:xudim] * cff2
    cff1 = 0.25 * (u[1:yvdim,:xvdim-2]+u[1:yvdim,1:xvdim-1]+u[:yvdim-1,:xvdim-2]+u[:yvdim-1,1:xvdim-1])
    cff2 = np.sqrt((v[1:yvdim,:xvdim-2]**2)+(cff1**2))
    bvstr[1:yvdim,:xvdim-2] = rdrg2 * v[1:yvdim,:xvdim-2] * cff2
    return bustr,bvstr

def ROMS_dz(nc,t=-999,x1=1,x2=-1,y1=1,y2=-1):
    H = nc.variables['h'][y1:y2,x1:x2]
    cs_w = nc.variables['Cs_w'][:]
    hc = nc.variables['hc'][:]
    sc_w = nc.variables['s_w'][:]
    Vtrans = nc.variables['Vtransform'][:]
    ydim, xdim = H.shape
    sdim = len(nc.dimensions['s_w'])
    if ( t < 0 ):
        zeta = 0.0
    else:
        zeta = nc.variables['zeta'][t,y1:y2,x1:x2]
    Z = np.zeros((sdim,(ydim),(xdim)))
    if (Vtrans == 1):
        for k in range(0,sdim):
            Z[k,:,:]   = (hc * sc_w[k] + (H-hc) * cs_w[k]) + (zeta * (1 + ((hc * sc_w[k] + (H-hc) * cs_w[k]) / H)))
    elif (Vtrans == 2):
        for k in range(0,sdim):
            Z[k,:,:]   = zeta + ((zeta + H) * ((hc * sc_w[k] + H * cs_w[k]) / (H+hc) ))
    else:
        print 'Unknown Vtrans - exit'
        sys.exit()
    dz=Z[1:,:,:]-Z[0:-1,:,:]
    return dz

def read_input():
    ''' read arguments from command line '''
    parser = argparse.ArgumentParser(description='Make curvilinear file for ROMS')
    ''' positional parameters '''
    parser.add_argument('yy',metavar='year', type=int, nargs=1,
                        help='Year')
    parser.add_argument('mm',metavar='month', type=int, nargs=1,
                        help='month')
    ''' options '''
    parser.add_argument('--hour',nargs='?',type=int,const=12,default=12,
                       help='hour (default value: 12)')
    args = parser.parse_args()
    return(args.yy[0],args.mm[0],args.hour)
