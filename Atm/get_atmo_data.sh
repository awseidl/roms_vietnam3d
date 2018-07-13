#!/bin/sh
## The contents of this file are based on the 2D operational script
## This script takes 3 input date arguments (YYYY, MM, DD) and downloads this from the specified server
## andrewws@met.no


## Server address, username and password
HOST='fog.oslo.dnmi.no'
USER=$(cat credentialsU)
PASSWD=$(cat credentialsP)

## Downloading atmospheric forcing file
echo "Downloading ECMWF atmospheric forcing file..."
echo "ECMWF file: vietnam_atmos_"$1$2$3"_12.nc"
pftp -i -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
binary
lcd $DIR_FIELDS
mget "vietnam_atmos_"$1$2$3"_12.nc"
quit
END_SCRIPT

TEMP_FILENAME="vietnam_atmo.nc"

mv vietnam_atmos_$1$2$3_12.nc $TEMP_FILENAME
echo "ECMWF data downloaded"

## Transform Atmo data
./Forcing_modification.bash $TEMP_FILENAME $4 $5
