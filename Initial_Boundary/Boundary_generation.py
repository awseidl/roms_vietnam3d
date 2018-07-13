## Properly produce Boundary files for 3D ROMS Vietnam East Sea. Input file must already be interpolated to curvilinear grid 
## This script takes 3 input arguments (input file (*_ocean.nc), output boundary condition file name, and ROMS working folder)
## andrewws@met.no

import netCDF4
import sys
import os

## Assigning input arguments
filein = sys.argv[1]
fileout = sys.argv[2]
outFolder = sys.argv[3]

print "Preparing Boundary conditions file for ROMS..."
#print "Input file: " +filein
print "Output file: " +fileout

## Determining and storing size of dimensions in input file
dimList = netCDF4.Dataset(filein)
xi_rho_num = int(dimList.dimensions['xi_rho'].size)
eta_rho_num = int(dimList.dimensions['eta_rho'].size)
xi_u_num = int(dimList.dimensions['xi_u'].size)
eta_v_num = int(dimList.dimensions['eta_v'].size)
s_rho_num = int(dimList.dimensions['s_rho'].size)

## Generating empty output netCDF with matching dimensions as input file
brync = netCDF4.Dataset(fileout, 'w', format='NETCDF4_CLASSIC')
s_rho = brync.createDimension('s_rho', s_rho_num)
eta_rho = brync.createDimension('eta_rho', eta_rho_num)
xi_rho = brync.createDimension('xi_rho', xi_rho_num)
xi_u = brync.createDimension('xi_u', xi_u_num)
eta_v = brync.createDimension('eta_v', eta_v_num)
ocean_time = brync.createDimension('ocean_time', None)

## Based of off the two variables listed below, extracts respective variable to respective boundary edge file, and renames them to something that ROMS can recognize. This file is then compressed (averaged) to eliminate the single value dimensions, before being concatenated to output boundary condition file
directions = ["North","South","East","West"]
input_vars = ["u","v","temp","salt","zeta","ubar","vbar"]
for bound in directions:
	print bound+ "ern boundary"
	if bound == "North" or bound == "South":
		dimensions = ["eta_rho","eta_v","eta_rho","eta_rho","eta_rho","eta_rho","eta_v"]
		if bound == "North":
			index = {'eta_rho': eta_rho_num-1,'eta_v': eta_v_num-1}
		else:
			index = {'eta_rho': 0,'eta_v': 0}
	else:
		dimensions = ["xi_u","xi_rho","xi_rho","xi_rho","xi_rho","xi_u","xi_rho"]
		if bound == "East":
			index = {'xi_rho': xi_rho_num-1,'xi_u': xi_u_num-1}
		else:
			index = {'xi_rho': 0,'xi_u': 0}
	for i,var_name in enumerate(input_vars,1):
		print var_name+ "_" +bound.lower()
		if i == 1:
			writeFlag = "-O"
		else:
			writeFlag = "-A" 
		boundSlice = "ncks " +writeFlag+ " -h -v " +var_name+ " -d " +dimensions[i-1]+ "," +str(index[dimensions[i-1]])+ " " +filein+ " " +fileout[:-3]+ "_" +bound.lower()+ ".nc" 
		boundRename = "ncrename -v " +var_name+ "," +var_name+ "_" +bound.lower()+ " " +fileout[:-3]+ "_" +bound.lower()+ ".nc"

		os.system(boundSlice)
		os.system(boundRename)
	dimsUsed = list(set(dimensions))
	## Averaging over single value dimensions
	for ii in dimsUsed:
		compression = "ncwa -O -a " +ii+ " " +fileout[:-3]+ "_" +bound.lower()+ ".nc " +fileout[:-3]+ "_" +bound.lower()+ ".nc"

		os.system(compression)
	## Concatenating file to output boundary condition file
	os.system("ncrcat -A " +fileout[:-3]+ "_" +bound.lower()+ ".nc " +fileout)
	os.system("rm " +fileout[:-3]+ "_" +bound.lower()+ ".nc")
	print "\n"

## Remove input file to reduce clutter and prevent cross contamination between simulations
os.system("rm " +filein)

## Moving boundary condition file into ROMS working folder
os.system("mv " +fileout+ " ../" +outFolder+ "/")
print "Boundary condition files ready in " +outFolder+ " folder"
