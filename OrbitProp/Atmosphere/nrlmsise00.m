function [out] = nrlmsise00(t,y,data)
% nrlmsise00 implements the NRLMSISE-00 atmospheric model. The output is
% compliant with the standard output.
% 
% Parameters:
%   data.nrlmsise00.day0 -> as the intial day of the year.
%
%   data.nrlmsise00.solar_activity to select different solar activity 
%   scenarios (options: 'Low', 'Medium', 'High', 'HighShort' or 'Custom'. 
%   If 'Custom' is selected then data.nrlmsise00.f107Average,
%   data.nrlmsise00.f107Daily and data.magneticIndex must be specified)
%
% Output:
%   out.Texo-> Exhospheric temperature in K.
%   out.Talt-> Altitude temperature in K.
%   out.rhoHe-> Helium particles per m3.
%   out.rhoO-> Atomic oxygen particles per m3.
%   out.rhoN2-> Diatomic nitrogen particles per m3.
%   out.rhoO2-> Diatomic oxygen particles per m3.
%   out.rhoAr-> Argon particles per m3.
%   out.rho-> HeTotal density in kg/m3.
%   out.rhoH-> Hydrogen particles per m3.
%   out.rhoN-> Atmoic Nitrogen particles per m3.
%   out.rhoOa-> Anomalous oxygen particles per m3.


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
dayOfYear = data.nrlmsise00.day0 + round(t/3600/24); %Day of the year
UTseconds = mod(t,3600*24);

%Solar activity
if strcmp(data.nrlmsise00.solar_activity,'Custom') %TODO allow solar input as function.
    f107Average = data.nrlmsise00.f107Average;
    f107Daily = data.nrlmsise00.f107Daily;
    magneticIndex = data.nrlmsise00.magneticIndex;
elseif strcmp(data.nrlmsise00.solar_activity,'Low')
    f107Average = 65;
    f107Daily = 65;
    magneticIndex = [0,0,0,0,0,0,0];
elseif strcmp(data.nrlmsise00.solar_activity,'Medium')
    f107Average = 140;
    f107Daily = 140;
    magneticIndex = [15,15,15,15,15,15,15];
elseif strcmp(data.nrlmsise00.solar_activity,'High')
    f107Average = 250;
    f107Daily = 250;
    magneticIndex = [45,45,45,45,45,45,45];
elseif strcmp(data.nrlmsise00.solar_activity,'HighShort')
    f107Average = 250;
    f107Daily = 300;
    magneticIndex = [240,240,240,240,240,240,240];
else 
    error('Error: Incorrect solar activity option')
end

%Latitude and longitude
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(y(1:3)',f,Req); %Compute taking assuming Earth as an ellipsoid
lat=lla(1);
lon=lla(2);
h=lla(3);

%Take into account rotation of the earth
lon = lon - t*we*180/pi;

%--- Compute the atmospheric properties ---%
[T rho] = atmosnrlmsise00(h, lat, lon, 0, dayOfYear, UTseconds, f107Average, f107Daily, magneticIndex,'None');
    
%--- Format density and temperature data ---%
out.Texo=T(:,1);
out.Talt=T(:,2);
out.rhoHe=rho(:,1);
out.rhoO=rho(:,2);
out.rhoN2=rho(:,3);
out.rhoO2=rho(:,4);
out.rhoAr=rho(:,5);
out.rho=rho(:,6);
out.rhoH=rho(:,7);
out.rhoN=rho(:,8);
out.rhoOa=rho(:,9);

end