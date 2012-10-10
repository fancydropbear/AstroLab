function [w] = WindHWM93(t,y,data)
% WindHWM93 implements the HWM93 wind models and gives the output in a 3x1
% wind vector in PCPF coordinates
%
% WARNING! This model needs compilation of Fortran code.
% Please make sure that the compiled version for your platform exists. If 
% not refer to the source code and compile.
%
% Parameters:
%   data.WindHWM93.day0 -> as the intial day of the year.
%
%   data.WindHWM93.solar_activity to select different solar activity 
%   scenarios (options: 'Low', 'Medium', 'High', 'HighShort' or 'Custom'. 
%   If 'Custom' is selected then data.WindHWM93.f107Average,
%   data.WindHWM93.f107Daily, data.WindHWM93.ap_daily and the
%   data.WindHWM93.ap_3hour must be specified).

%--- Copyright notice ---%
% Copyright 2012 Cranfield University
% Written by Josep Virgili
%
% This file is part of the AstroLab toolkit
%
% AstroLab is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% AstroLab is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with AstroLab.  If not, see <http://www.gnu.org/licenses/>.

%---- CODE ----%

%--- Prepeare input data --%
%Date and time
iyd = data.WindHWM93.day0 + round(t/3600/24); %Day of the year
sec = mod(t,3600*24);
slt = sec/3600; %Local aparent solar time TODO

%Latitude and longitude
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(y(1:3)',f,Req); %Compute taking assuming Earth as an ellipsoid
lat=lla(1);
lon=lla(2);
alt=lla(3)/1e3;

%Solar activity
if strcmp(data.WindHWM93.solar_activity,'Custom') %TODO allow solar input as function.
    f107a = data.WindHWM93.f107Average;
    f107 = data.WindHWM93.f107Daily;
    ap_daily = data.WindHWM93.ap_daily;
    ap_3hour = data.WindHWM93.ap_3hour;
elseif strcmp(data.nrlmsise00.solar_activity,'Low')
    f107a = 65;
    f107 = 65;
    ap_daily = 0;
    ap_3hour = 0;
elseif strcmp(data.nrlmsise00.solar_activity,'Medium')
    f107a = 140;
    f107 = 140;
    ap_daily = [15,15,15,15,15,15,15];
    ap_3hour = 15;
elseif strcmp(data.nrlmsise00.solar_activity,'High')
    f107a = 250;
    f107 = 250;
    ap_daily = 45;
    ap_3hour = 45;
elseif strcmp(data.nrlmsise00.solar_activity,'HighShort')
    f107a = 250;
    f107 = 300;
    ap_daily = 240;
    ap_3hour = 240;
else 
    error('Error: Incorrect solar activity option')
end

%--- Compute wind ---%
wind = hwm93_raw(iyd,sec,alt,lat,lon,slt,f107a,f107,ap_daily,ap_3hour); 

%--- Convert to ECEF frame ---%

%Take into account the rotation of the earth
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
theta = t*we;
Rot = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];

%Coonvert to ECEF frame
w = Rot\[-wind(2)*sind(lat);wind(1);wind(2)*cosd(lat)];


end