#!/bin/bash
## Properly setup CMEMS initialisation file for 3D ROMS and call functions to make boundary file.
## This script takes 6 arguments (YYYY, MM, DD, name of ocean state initialisation file, boundary condition file name, and ROMS working folder 
## andrewws@met.no

echo "Preparing Initial ocean state file for ROMS..."
FILEIN=vietnam3d_forecast.nc
FILEOUT_TEMP=vietnam3d_forecast_with2D.nc
#echo $FILEIN
#echo $FILEOUT_TEMP

#FILEORIGINAL="${FILEIN::-3}_Original.nc"
#cp $FILEIN $FILEORIGINAL
#echo Original copied

## Renaming variables to names that ROMS will recognize
ncrename -v thetao,temp $FILEIN
ncrename -v so,salt $FILEIN
ncrename -v vo,v $FILEIN
ncrename -v uo,u $FILEIN
ncrename -v zos,zeta $FILEIN
echo "Variables renamed"

## Determine thickness of layers
ncap2 -O -s "thckns=depth(1:depth.size()-1)-depth(0:depth.size()-2)" $FILEIN temp.nc
ncap2 -O -h -s "thckns(0)=0" temp.nc temp.nc

## Calculate 2D momentums
echo "Calculating 2D momentums...u"
# u
ncwa -O -a depth -w thckns -v u temp.nc temp2.nc
ncrename -v u,ubar temp2.nc
ncks -A -h -v ubar temp2.nc $FILEIN
echo "done...v"
## v
ncwa -O -a depth -w thckns -v v temp.nc temp2.nc
ncrename -v v,vbar temp2.nc
ncks -A -h -v vbar temp2.nc $FILEIN
echo "done"
rm temp.nc
rm temp2.nc

## Modify variable attributes
for v in temp salt u v ubar vbar; do 
	ncatted -O -h -a standard_name,$v,d,, $FILEIN; 
done
ncatted -O -h -a long_name,ubar,m,c,"vertically integrated u-momentum component" $FILEIN
ncatted -O -h -a long_name,vbar,m,c,"vertically integrated v-momentum component" $FILEIN
ncatted -O -h -a long_name,temp,m,c,"sea water potential temperature" $FILEIN
ncatted -O -h -a long_name,salt,m,c,"sea water salinity" $FILEIN
ncatted -O -h -a long_name,u,m,c,"u-momentum component" $FILEIN
ncatted -O -h -a long_name,v,m,c,"v-momentum component" $FILEIN
ncatted -O -h -a long_name,zeta,m,c,"free surface" $FILEIN
echo "Variable attributes modified"

## Set reference time to days since UNIX epoch
cdo setreftime,1970-1-1,00:00:00,days $FILEIN $FILEOUT_TEMP
echo "Ref Time set to: Days since 1970-01-01 00:00"

## Remove raw input file to reduce clutter and prevent cross contamination between simulations
rm $FILEIN


######### Grid interpolation ############
FILEOUT=vietnam3d_ocean.nc
#echo $FILEOUT_TEMP
#echo $FILEOUT

## Python script provided by nilsmk@met.no for interpolating CMEMS data to ROMS curvilinear grid
python make_vietnam3d_opera.py $FILEOUT_TEMP

## Copying ocean state initialisation file (*_ini.nc) into ROMS folder
cp Output/"vietnam3d_01.nc" ../$6/$4
echo 'Initial ocean state file ready in "'$6'" folder'

## Moving into Output folder for working with interpolated files
cd Output
## Number of daily files
NFILES=$(ls -1 | wc -l)

## concatenating netCDF files into one
ncrcat -n $NFILES,2,1 vietnam3d_01.nc $FILEOUT 

## Moving FILEOUT up one folder level
mv $FILEOUT ../${FILEOUT}

## Removing individual netCDF files in Output folder to reduce clutter and prevent cross-contamination between simulations
rm vietnam3d_*

## Moving out of Output folder
cd ..

## Remove modified input file to reduce clutter and prevent cross contamination between simulations
rm $FILEOUT_TEMP

## Produce boundary conditions (*_bry.nc) from ocean file (*_ocean.nc)
python Boundary_generation.py $FILEOUT $5 $6
#./Boundary_production.bash $FILEOUT $5 $6



