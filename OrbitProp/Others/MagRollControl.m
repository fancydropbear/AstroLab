function [A,T,E] = MagRollControl(t,y,data)
% MagRollControls is a model that reduces to 0 the roll by aplying a torque
% proportional to the roll angle. The torque is assumed to be provided by
% magnetorquers
%
% The Earth magnetic field is modeled using the IGRF-11 and therfore the
% model requires the igrf11magm MATLAB function.
%
% Parameters:
%   data.MagRollControl.A -> Actuation level in Am2
%
%   data.MagRollControl.Y -> Year (include decimal point for increased precision)
%
%   data.MagRollControl.k -> Proportional constant. T = roll * k
%
%   data.MagRollControl.wind_model -> Function handle specifying the wind
%   model. The output must be a 3x1 vector containg the wind in PCPF
%   coordinates.
%

%--- Copyright notice ---%
% Copyright 2012 Cranfield University
% Written by Josep Virgili
%
% This file is part of the AstroLab
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

%--- CODE ---%

%--- Get magnetic field ---%
%Latitude and longitude
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(y(1:3)',f,Req); %Compute taking assuming Earth as an ellipsoid
lat=lla(1);
lon=lla(2);
h=lla(3);
%Take into account rotation of the earth
lon = mod(lon - t*we*180/pi,360);
if lon>180
    lon = lon-360;
elseif lon<-180
    lon = lon+360;
end

%Magnetic field in T
mfield = igrf11magm(h, lat, lon, data.MagRollControl.Y)/1e9;

%--- Max torque ---%
Tmax = abs(cross(mfield,data.MagRollControl.A));

%--- Attitude in flow reference frame ---%
%Flow is velocity - atmosphere rotatio - wind
wind = data.MagRollControl.wind_model(t,y,data); %Wind
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
corot = cross([0;0;we],y(1:3)); %Atmosphere rotation
flow = y(4:6) - corot - wind;
%Flow reference frame
xa = flow/norm(flow); %Velocity (roll axis)
zaux = -y(1:3)/norm(y(1:3)); %Nadir pointing
ya = cross(zaux,xa); % Velocity anti-paralell (pitch axis)
za = cross(xa,ya); % Near nadir (yaw axis)

%Rotation matrix from ECEF to flow ref
Rfe=[xa';ya';za'];

%Rot from ECEF to body
Rbe=quat2dcm(y(7:10)');

%Rot from flow to body
Rbf = Rbe*inv(Rfe);
%Get angles
[r1, r2, r3] = dcm2angle(Rbf,'XYZ');

%---- Torque ---%
T=[-r1*data.MagRollControl.k;0;0];

%Check if desired Torque is not bigger than the max
if T(1)>Tmax(1)
    %Actuator saturated
    if data.MagRollControl.verb; disp('X actuator saturating'); end;
    T(1)=sign(T(1))*Tmax(1);
end

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);


end