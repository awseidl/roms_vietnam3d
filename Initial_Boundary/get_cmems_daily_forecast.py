## This script takes 4 input arguments (YYYY, MM, DD, number of days in simulation)
## It then automatically determines the appropriate ending date for download and downloads CMEMS data according to specifications below
## andrewws@met.no

from datetime import datetime, timedelta
import sys
import os

## Assigning input arguments
startY = sys.argv[1]
startM = sys.argv[2]
startD = sys.argv[3]
simDays = int(sys.argv[4])

## Start and end dates as datetimes 
startDate = datetime(int(startY),int(startM),int(startD))
endDate = startDate + timedelta(simDays)
endY = str(endDate.year)
endM = str("%02d" % endDate.month)
endD = str("%02d" % endDate.day)

key = open("credentials","r")
credentials = key.read().split("\n")

## Command to call in terminal
command = "python motu-client-python/motu-client.py " +\
"--user '" + credentials[0] + "' " +\
"--pwd '" + credentials[1] + "' " +\
"--motu http://nrt.cmems-du.eu/motu-web/Motu " +\
"--service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS " +\
"--product-id global-analysis-forecast-phy-001-024 " +\
"--longitude-min 98 " +\
"--longitude-max 122 " +\
"--latitude-min 3 " +\
"--latitude-max 24 " +\
"--date-min "+startY+"-"+startM+"-"+startD+" 12:00:00 " +\
"--date-max "+endY+"-"+endM+"-"+endD+" 12:00:00 " +\
"--depth-min 0.493 " +\
"--depth-max 5727.918000000001 " +\
"--variable thetao --variable so --variable zos --variable uo --variable vo " +\
"--out-dir . " +\
"--out-name vietnam3d_forecast.nc"

print "Downloading CMEMS ocean state file..."

os.system(command)

## Original contents of get_cmems_daily_forecast text file (and a description of what they represent)
#/usr/bin/python \
#motu-client-python/motu-client.py \ 				# filepath to motu python file (from http://marine.copernicus.eu/faq/what-are-the-motu-and-python-requirements/)
#--user '' \	 						# username
#--pwd ''	 \ 					# plaintext password
#--motu http://nrt.cmems-du.eu/motu-web/Motu \			# MOTU server
#--service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS \ 	# model used
#--product-id global-analysis-forecast-phy-001-024 \		# also model used
#--longitude-min 98 \ 						# westernmost point
#--longitude-max 122 \						# easternmost point
#--latitude-min 3 \ 						# southernmost point
#--latitude-max 24 \ 						# northernmost point
#--date-min "$1-$2-$3 12:00:00" \ 				# starting time
#--date-max "$4-$5-$6 12:00:00" \				# ending time
#--depth-min 0.493 \						# highest depth
#--depth-max 5727.918000000001 \ 				# deepest depth
#--variable thetao --variable so --variable zos --variable uo --variable vo \ # potential temperature [C], salinity [psu], sea surface height [m], eastward velocity [m/s], northward velocity [m/s]
#--out-dir . \							# directory to save output
#--out-name vietnam3d_forecast.nc" 				# output file name
