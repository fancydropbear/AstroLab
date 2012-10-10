clear;

iyd = 172;
sec = 29000;
alt = 400;
lat = 60;
lon = -70;
slt = 16;
f107a = 150;
f107 = 150;
ap_daily = 4;
ap_3hour = 4;

my_winds = hwm93(iyd,sec,alt,lat,lon,slt,f107a,f107,ap_daily,ap_3hour)

% compare with sample_output_hwm93.txt