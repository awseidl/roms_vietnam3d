# This file automatically changes text in the ".in" (input) file! 
# Please check this file carefully if you find things are unwantingly changing in the .in file and you don't know why
# andrewws@met.no

import netCDF4
import sys
from datetime import date
import os

inputFile = sys.argv[1]
atmFile = sys.argv[2]
initialFile = sys.argv[3]
bryFile = sys.argv[4]
NCFolder = sys.argv[5]
rootFolder = sys.argv[6]
simDays = int(sys.argv[7])

ini = netCDF4.Dataset(NCFolder+ "/" +initialFile)
startTime = ini.variables['ocean_time'][0]
startDate = netCDF4.num2date(startTime,ini.variables['ocean_time'].units)

# Updating start date
dateReplace = "      DSTART =  " + str(startTime) + "d0                  ! " + str(startDate) + "\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if line.find("DSTART =") != -1:
            f.writelines(dateReplace)
        else:
            f.writelines(line)

# Updating Initial ocean state netcdf path
iniReplace = "   ININAME == " +rootFolder+ "/" +NCFolder+ "/" +initialFile+ "        ! NLM initial conditions\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if line.find("ININAME ==") != -1:
            f.writelines(iniReplace)
        else:
            f.writelines(line)

# Updating Boundary condition netcdf path
bryReplace = "   BRYNAME == " +rootFolder+ "/" +NCFolder+ "/" +bryFile+ "\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("BRYNAME ==") != -1) & (line[0] != "!"):
            f.writelines(bryReplace)
        else:
            f.writelines(line)

# Updating Atmospheric forcing netcdf path
atmReplace = "   FRCNAME == " +rootFolder+ "/" +NCFolder+ "/" +atmFile+ "\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("FRCNAME ==") != -1) & (line[0] != "!"):
            f.writelines(atmReplace)
        else:
            f.writelines(line)

# Updating date auto used
usedReplace = "!       -Last auto use:  " +str(date.today())+ "                                           !\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("-Last auto use:") != -1):
            f.writelines(usedReplace)
        else:
            f.writelines(line)

# Looking for VARNAME in correct place
varnameReplace = "   VARNAME = " +rootFolder+ "/trunk/ROMS/External/varinfo.dat\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("VARNAME =") != -1):
            f.writelines(varnameReplace)
        else:
            f.writelines(line)

# Looking for grid in correct place
grdReplace = "   GRDNAME == " +rootFolder+ "/Initial_Boundary/east_sea_grd-smooth.nc   ! Grid\n"
with open(inputFile,'r') as f:
    get_all=f.readlines()

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("GRDNAME ==") != -1):
            f.writelines(grdReplace)
        else:
            f.writelines(line)

# Setting timesteps according to user specified simulation length
with open(inputFile,'r') as f:
    get_all=f.readlines()
    for i,line in enumerate(get_all,1):  
        if (line.find("DT ==") != -1):
            DT = int(line[line.rfind('=')+1:])
 
NTimesReplace = "      NTIMES == " +str(60/DT*1440*simDays) + "   ! " +str(simDays)+ " days at DT=" +str(DT)+ " seconds\n"

with open(inputFile,'w') as f:
    for i,line in enumerate(get_all,1):  
        if (line.find("NTIMES ==") != -1):
            f.writelines(NTimesReplace)
        else:
            f.writelines(line)

print "Changes made to " +inputFile
