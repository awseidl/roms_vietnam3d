#!/bin/bash
## Properly setup ECMWF forcing file for 3D ROMS
## This script take 3 input arguments (input file, output file, ROMS working folder)
## andrewws@met.no

FILEIN=$1
FILEOUT=$2
echo "Preparing Atmospheric forcing file for ROMS..."
echo "Output file: " $FILEOUT

#FILEORIGINAL="${FILEIN::-3}_Original.nc"
#cp $FILEIN $FILEORIGINAL
#echo Original copied

## There are many fields present in file that we will not be using. They are removed to reduce file size and speed up upcoming netCDF operations
ncks -O -h -x -v fog_area_fraction,air_temperature_pl,relative_humidity_pl,y_wind_pl,x_wind_pl,geopotential_pl,surface_geopotential,low_type_cloud_area_fraction,surface_air_pressure,medium_type_cloud_area_fraction,land_area_fraction,lwe_thickness_of_surface_snow_amount,sea_surface_temperature,high_type_cloud_area_fraction,soil_temperature,lwe_thickness_of_convective_precipitation_amount,surface_downwelling_shortwave_flux_in_air_acc,lwe_thickness_of_stratiform_precipitation_amount,specific_convective_available_potential_energy,surface_upward_latent_heat_flux_acc,wind_speed_of_gust,surface_upward_sensible_heat_flux_acc,ap,b,depth_between_layers,forecast_reference_time,hybrid,p0,pressure,surface_air_pressure $FILEIN $FILEIN

## ECMWF forcing file comes flipped in North-South direction (Min latitude 40N, max latitude 0). This command flips all fields and the dimension correctly
ncpdq -O -a -latitude $FILEIN $FILEIN
echo "Latitude flipped"

## Renaming variables of interest to names that ROMS will recognize
ncrename -d longitude,lon -v longitude,lon $FILEIN
ncrename -d latitude,lat -v latitude,lat $FILEIN
ncrename -v y_wind_10m,Vwind $FILEIN
ncrename -v x_wind_10m,Uwind $FILEIN
ncrename -v cloud_area_fraction,cloud $FILEIN
ncrename -v air_temperature_2m,Tair $FILEIN
ncrename -v dew_point_temperature_2m,Tdair $FILEIN
ncrename -v air_pressure_at_sea_level,Pair $FILEIN
echo "Variables renamed"

## Modify time attributes
ncatted -O -h -a long_name,time,m,c,time $FILEIN
ncatted -O -h -a standard_name,time,c,c,time $FILEIN
echo "Time attributes modified"

## Convert values: pressures to hPa
ncap2 -O -s "Pair=Pair/100" $FILEIN $FILEIN
ncatted -O -h -a units,Pair,m,c,"hPa" $FILEIN
ncatted -O -h -a long_name,Pair,m,c,"Mean sea level pressure" $FILEIN
echo "Pressure converted to hPa"

## Convert values: temps to Celcius
ncap2 -O -s "Tair=Tair-273.16" $FILEIN $FILEIN
ncatted -O -h -a units,Tair,m,c,"C" $FILEIN
ncatted -O -h -a long_name,Tair,m,c,"2 metre temperature" $FILEIN

ncap2 -O -s "Tdair=Tdair-273.16" $FILEIN $FILEIN
ncatted -O -h -a units,Tdair,m,c,"C" $FILEIN
ncatted -O -h -a long_name,Tdair,m,c,"2 metre dewpoint temperature" $FILEIN
echo "Temperatures converted to Celsius"

## Calculate relative humidity according to Tetens" formula (cf bulk_flux.F).
ncap2 -O -s "Qair = 100*exp(17.502*( (Tdair/(240.97+Tdair)) - (Tair/(240.97+Tair)) ))" $FILEIN $FILEIN
ncatted -O -h -a units,Qair,m,c,percentage $FILEIN
ncatted -O -h -a long_name,Qair,c,c,"2 metre relative humidity" $FILEIN
ncatted -O -h -a standard_name,Qair,m,c,"air_relative_humidity" $FILEIN
echo "RH calculated"

## Set reference time to days since UNIX epoch
cdo setreftime,1970-1-1,00:00:00,days $FILEIN $FILEOUT
echo "Ref Time set to: Days since 1970-01-01 00:00"

## Modify variable attributes
for v in Vwind Uwind cloud Tair Tdair Pair Qair; do
	ncatted -O -h -a coordinates,$v,c,c,"lon lat" $FILEOUT
done
ncwa -O -a surface $FILEOUT $FILEOUT
echo "Variable attributes modified"


## Remove raw input forcing file to reduce clutter and prevent cross contamination between simulations
rm $FILEIN

## For some reason, 26th (25th with 0 indexing) time is always composed of entirely missing data. This ncks command removes it
ncks -O -d time,0,24 -d time,26, $FILEOUT $FILEOUT

mv $FILEOUT ../$3/
echo 'Atmospheric forcing file ready in "'$3'" folder\n'