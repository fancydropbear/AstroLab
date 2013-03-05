function [A,T,E] = MagDamping(t,y,data)
% MagDamping is a model that damps the unwanted angular rates by using 
% magnetic torquer rods.
%
% The Earth magnetic field is modeled using the IGRF-11 and therfore the
% model requires the igrf11magm MATLAB function.
%
% Parameters:
%   data.MagDamping.c -> Torque per rad/s. Damping ratio
%
%   data.MagDamping.A -> Actuation level vector in Am2
%
%   data.MagDamping.Y -> Year (include decimal point for increased precision)
%
%   data.MagDamping.w -> Body rate bias rad/s
%

%--- Copyright notice ---%
% Copyright 2012-2013 Cranfield University
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
lon = mod(lon - t*we*180/pi,365);
if lon>180
    lon = lon-360;
elseif lon<-180
    lon = lon+360;
end

%Magnetic field in T
mfield_ECEF = igrf11magm(h, lat, lon, data.MagDamping.Y)/1e9;
%Change magnetic field to body axes.
DCM = quat2dcm([y(10),y(7:9)']);
mfield = DCM*mfield_ECEF';

%--- Desired torque ---%
Tdes = -(y(11:13)-data.MagDamping.w')*data.MagDamping.c;

%--- Check if torque is possible ---%

%Preliminary actuation levels m_p
m_p=(cross(mfield,Tdes)/norm(mfield)^2)';
%Check different scalar values
k_max=(data.MagDamping.A-m_p)./mfield'; %Scalar value for maximum actuation level on one axis
k_min=(-data.MagDamping.A-m_p)./mfield'; %Scalar value for minus maximum actuation level on one axis
k_0=-m_p./mfield'; %Scalar value for 0 actuation level on one axis

%Scan through the different scalar values
err=1e-3/100;
for k=[k_max,k_min,k_0]
    m=m_p+k*mfield';
    
    if abs(m(1))>data.MagDamping.A(1)*(1+err)
        %Solution is not real
        continue
    end
    if abs(m(2))>data.MagDamping.A(2)*(1+err)
        %Solution is not real
        continue
    end
    if abs(m(3))>data.MagDamping.A(3)*(1+err)
        %Solution is not real
        continue
    end
    
    %Solution is between bounds
    T=Tdes;
    
    %Acceleration (this model don't produce linear acceleration)
    A=[0;0;0];
    
    %Extra state variables (this model doesn't need extra state variables)
    E=zeros(1,length(y)-13);
    
    return
end

%Torque is not achievable, actuators saturate, get the nearest one.
if data.MagDamping.verb; disp('Magnetic torquers saturating'); end;

%Scan through the different scalar values
T_pos=[];
T_diff=[];
for k=[k_max,k_min,k_0]
    m=m_p+k*mfield';
    
    if abs(m(1))>data.MagDamping.A(1)
        %Adjust to maximum
        m(1)=sign(m(1))*data.MagDamping.A(1);
    end
    if abs(m(2))>data.MagDamping.A(2)
        %Adjust to maximum
        m(2)=sign(m(2))*data.MagDamping.A(2);
    end
    if abs(m(3))>data.MagDamping.A(3)
        %Adjust to maximum
        m(3)=sign(m(3))*data.MagDamping.A(3);
    end
    
    %Compute possible torque
    T_pos(end+1,1:3)=cross(m,mfield);
    %Compute difference to desired torque
    T_diff(end+1)=norm(T_pos(end,:)'-Tdes);
    
end

%Get the configuration with minimum difference
[~,I]=min(T_diff);
T=T_pos(I,:)';

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end