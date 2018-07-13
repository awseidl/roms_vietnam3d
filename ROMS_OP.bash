#!/bin/bash
###################### Information about this package ##################################
## Most of these scripts were designed using NCO version 4.5.4, and CDO version 1.7.0 ##
## Please contact andrewws@met.no or lrh@met.no for more information regarding use    ##
########################################################################################
  
##### Inputs #####
## Simulation duration (days). ECMWF limitations mean this can be a maximum of 10! 
SIM_DAYS=5

## Filenames
ROMS_ATM_FILENAME="vietnam3d_atm.nc"
ROMS_INI_FILENAME="vietnam3d_ini.nc"
ROMS_BRY_FILENAME="vietnam3d_bry.nc"
ROMS_INPUT_FOLDER="ROMS_Input-Output"
ROMS_INPUT_FILE="vietnam3d.in"
ROMS_LOG_FILE="vietnam3d.log"


###################### Automated things below ##########################################
## Automatically set start date to the previous day (critical for ECMWF download) 
SYEAR=`date '+%Y' -d "-1 day"`
SMONTH=`date '+%m' -d "-1 day"`
SDAY=`date '+%d' -d "-1 day"`

## Finding location of this folder on system
ROOT_FOLDER=${PWD}

## Reviewing 
HOLD_ON=true
while $HOLD_ON; do
echo "Simulation start time: "$SYEAR"-"$SMONTH"-"$SDAY" 12:00"
echo "Simulation end time:   "`date '+%Y-%m-%d' -d "+$(( SIM_DAYS - 1 )) day"`" 12:00"
echo "Simulation duration:            "$SIM_DAYS" days"
echo "Run this simulation?" 
echo "(N/n: change simulation duration; All other input: starts model run)"
read answer

## Checking for valid inputs
if [ "$answer" != "${answer#[Nn]}" ]; then
    echo "Please enter desired simulation duration (whole days, 1 to 10):"
    read SIM_DAYS
    INVALID_NUM=true
    while $INVALID_NUM; do
        if ! [[ "$SIM_DAYS" =~ ^[\-0-9]+$ ]]; then
            echo "Integers required..."
            read SIM_DAYS
        elif [[ $SIM_DAYS -lt 0 ]]; then
            echo "Greater than 0..."
            read SIM_DAYS
        elif [[ $SIM_DAYS -gt 10 ]]; then
            echo "Maximum of 10..."
            read SIM_DAYS
        else
            INVALID_NUM=false
        fi
    done
else
    echo "Engage!"
    HOLD_ON=false
fi
done

## Download Atmo forcings and transform as neccesary
cd ${ROOT_FOLDER}/Atm
./get_atmo_data.sh $SYEAR $SMONTH $SDAY $ROMS_ATM_FILENAME $ROMS_INPUT_FOLDER

## Download Ocean state and forecast
cd ${ROOT_FOLDER}/Initial_Boundary
python get_cmems_daily_forecast.py $SYEAR $SMONTH $SDAY $SIM_DAYS

## Preparing Initial ocean state and Boundary conditions as neccessary
./Initial_modification.bash $SYEAR $SMONTH $SDAY $ROMS_INI_FILENAME $ROMS_BRY_FILENAME $ROMS_INPUT_FOLDER

cd ${ROOT_FOLDER}
## Build ROMS executable (oceanM)
echo "Building ROMS..."
./build_vietnam3d.bash -j 2

## Prepare input (*.in) file for ROMS
python prepareInputFile.py $ROMS_INPUT_FILE $ROMS_ATM_FILENAME $ROMS_INI_FILENAME $ROMS_BRY_FILENAME $ROMS_INPUT_FOLDER $ROOT_FOLDER $SIM_DAYS

## Prepare the viewing of empty log file
rm $ROMS_LOG_FILE
echo > $ROMS_LOG_FILE
gnome-terminal -e "tail -f "$ROMS_LOG_FILE

## Start ROMS
echo "Now running ROMS..."
mpirun -n 2 oceanM $ROMS_INPUT_FILE >> $ROMS_LOG_FILE

## Visual cue
for i in {0..5}; do echo "##      ##"; done

## Copying completed output to separate folder
echo "Simulation complete"
cp $ROMS_LOG_FILE Complete_output/${ROMS_LOG_FILE}
cd ${ROMS_INPUT_FOLDER}
cp vietnam3d_his.nc ../Complete_output/vietnam3d_his.nc
cp vietnam3d_avg.nc ../Complete_output/vietnam3d_avg.nc
cp vietnam3d_sta.nc ../Complete_output/vietnam3d_sta.nc
echo 'Output netCDFs can be found in "Complete_output" folder'


