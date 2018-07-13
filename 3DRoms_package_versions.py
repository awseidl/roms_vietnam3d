#!/usr/bin/python
##

import os

packages = ['gfortran','netcdf_c','netcdf_fortran','nco','cdo','flex','m4','openmpi','pip','curl','ftp']

versions_req = {'gfortran': "4:5.3.1", \
		'netcdf_c': "4.6.1", \
		'netcdf_fortran': "4.4.4", \
		'nco': "4.5.4", \
		'cdo': "1.7.0", \
		'flex': "2.6.0", \
		'm4': "1.4.17", \
		'openmpi': "3.0.1", \
		'pip': "8.1.1", \
		'curl': "7.47.0", \
		'ftp': "0.17"}

## gfortran
try:
	gfortran = os.popen("dpkg -s gfortran | grep -i version").read()
	gfortran = gfortran.split(" ")[1]
	gfortran = gfortran.split("\n")[0]
	gfortran = gfortran.split("+")[0]
	gfortran = gfortran.split("-")[0]
except IndexError:
	gfortran = "Missing"

## NetCDF-C
try:
	netcdf_c = os.popen("nc-config --version").read()
	netcdf_c = netcdf_c.split(" ")[1]
	netcdf_c = netcdf_c.split("\n")[0]
except IndexError:
	netcdf_c = "Missing"

## NetCDF-Fortran
try:
	netcdf_f = os.popen("nf-config --version").read()
	netcdf_f = netcdf_f.split(" ")[1]
	netcdf_f = netcdf_f.split("\n")[0]
except IndexError:
	netcdf_f = "Missing"
	
## NCO
try:
	nco = os.popen("dpkg -s nco | grep -i version").read()
	nco = nco.split(" ")[1]
	nco = nco.split("\n")[0]
	nco = nco.split("+")[0]
	nco = nco.split("-")[0]
except IndexError:
	nco = "Missing"

## CDO
try:
	cdo = os.popen("dpkg -s cdo | grep -i version").read()
	cdo = cdo.split(" ")[1]
	cdo = cdo.split("\n")[0]
	cdo = cdo.split("+")[0]
	cdo = cdo.split("-")[0]
except IndexError:
	cdo = "Missing"

## flex
try:
	flex = os.popen("dpkg -s flex | grep -i version").read()
	flex = flex.split(" ")[1]
	flex = flex.split("\n")[0]
	flex = flex.split("+")[0]
	flex = flex.split("-")[0]
except IndexError:
	flex = "Missing"

## m4
try:
	m4 = os.popen("dpkg -s m4 | grep -i version").read()
	m4 = m4.split(" ")[1][:-1]
	m4 = m4.split("\n")[0]
	m4 = m4.split("+")[0]
	m4 = m4.split("-")[0]
except IndexError:
	m4 = "Missing"

## OpenMPI
try:
	openmpi = os.popen("mpirun --version").read()
	openmpi = openmpi.split(" ")[3]
	openmpi = openmpi.split("\n")[0]
except IndexError:
	openmpi = "Missing"

## pip
try:
	pip = os.popen("pip --version").read()
	pip = pip.split(" ")[1]
except IndexError:
	pip = "Missing"

## curl
try:
	curl = os.popen("dpkg -s curl | grep -i version").read()
	curl = curl.split(" ")[1]
	curl = curl.split("\n")[0]
	curl = curl.split("+")[0]
	curl = curl.split("-")[0]
except IndexError:
	curl = "Missing"

## ftp
try:
	ftp = os.popen("dpkg -s ftp | grep -i version").read()
	ftp = ftp.split(" ")[1]
	ftp = ftp.split("\n")[0]
	ftp = ftp.split("+")[0]
	ftp = ftp.split("-")[0]
except IndexError:
	ftp = "Missing"


versions_install = {'gfortran': gfortran, \
		'netcdf_c': netcdf_c, \
		'netcdf_fortran': netcdf_f, \
		'nco': nco, \
		'cdo': cdo, \
		'flex': flex, \
		'm4': m4, \
		'openmpi': openmpi, \
		'pip': pip, \
		'curl': curl, \
		'ftp': ftp}

#sys.exit()
print "Package/Software Name\t\tVersion:"
print "\t\t\tRecommended\tInstalled"
for i in packages:
	if len(i) <= 7:
		spacer = "\t\t\t"
	else:
		spacer = "\t\t"
	print i+ spacer + versions_req[i] + "\t\t" + versions_install[i]
