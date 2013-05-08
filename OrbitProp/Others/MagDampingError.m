function [A,T,E] = MagDampingError(t,y,data)
% MagDampingError is a explicit model that damps the unwanted angular rates by using 
% magnetic torquer rods.It is an extensional version of MagDamping.m file.
% Error/Bias/Noise source from all instruments are consideraded in this
% case.
%
% The Earth magnetic field is modeled using the IGRF-11 and therfore the
% model requires the igrf11magm MATLAB function.
%
% Parameters:
%   data.MagDampingError.c -> Torque per rad/s. Damping ratio
%
%   data.MagDampingError.A -> Actuation level vector in Am2
%
%   data.MagDampingError.Y -> Year (include decimal point for increased precision)

%   data.MagDampingError.w -> Body rate bias rad/s

%   data.MagDampingError.MTB -> Actuation Level Bias Vector in Am2

%   data.MagDampingError.MMR -> Magnetometre Resolution in Teslas 

%   data.MagDampingError.MTR -> Magnetorquer Resolution in Am2 

%   data.MagDampingError.GE -> Gyroscope Reading One Sigma Uncertaincy Vector in rad/s

%   data.MagDampingError.SF ->Gyroscope Sampling/Error Frequency in Hz

%--- Copyright notice ---%
% Copyright 2012-2013 Cranfield University
% Written by Josep Virgili and Daniel Zhou Hao
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
mfield_ECEF = igrf11magm(h, lat, lon, data.MagDampingError.Y)/1e9;
%Change magnetic field to body axes.
DCM = quat2dcm([y(10),y(7:9)']);
mfield = DCM*mfield_ECEF';
mfield_round = round(mfield./data.MagDampingError.MMR).*data.MagDampingError.MMR;

%--- Desired torque ---%
% angular_rates = normrnd (y(11:13), data.MagDampingError.GE);
% This random angular_rate induced stiff problem to the integrator, make
% the progress slow

% Bias Appoarch 
angular_rates = y(11:13) + data.MagDampingError.GE .* sin(2*pi*data.MagDampingError.SF*t); % 25 HZ sampling frq.
Tdes = -(angular_rates - data.MagDampingError.w')*data.MagDampingError.c;

%display time for debugging
% t/60 ;

%--- Check if torque is possible ---%

%Preliminary actuation levels m_p
m_p=(cross(mfield_round,Tdes)/norm(mfield_round)^2)';

%Check different scalar values
k_max=(data.MagDampingError.A-m_p)./mfield_round'; %Scalar value for maximum actuation level on one axis
k_min=(-data.MagDampingError.A-m_p)./mfield_round'; %Scalar value for minus maximum actuation level on one axis
k_0=-m_p./mfield_round'; %Scalar value for 0 actuation level on one axis

%Scan through the different scalar values
err=1e-3/100;
for k=[k_max,k_min,k_0]
    m=m_p+k*mfield_round';
  
    if abs(m(1))>data.MagDampingError.A(1)*(1+err)
        %Solution is not real
        continue
    end
    if abs(m(2))>data.MagDampingError.A(2)*(1+err)
        %Solution is not real
        continue
    end
    if abs(m(3))>data.MagDampingError.A(3)*(1+err)
        %Solution is not real
        continue
    end
    
    m = round(m(:)./data.MagDampingError.MTR(:)).*data.MagDampingError.MTR(:)+data.MagDampingError.MTB(:); 
    
    %Solution is between bounds
    T=cross(m(:),mfield(:));
    
    %Acceleration (this model don't produce linear acceleration)
    A=[0;0;0];
    
    %Extra state variables (this model doesn't need extra state variables)
    E=zeros(1,length(y)-13);
    
    return
end

%Torque is not achievable, actuators saturate, get the nearest one.
if data.MagDampingError.verb; disp('Magnetic torquers saturating'); end;

%Scan through the different scalar values
T_pos=[];
m_pos=[];
T_diff=[];
for k=[k_max,k_min,k_0]
    m=m_p+k*mfield_round';
    
    if abs(m(1))>data.MagDampingError.A(1)
        %Adjust to maximum
        m(1)=sign(m(1))*data.MagDampingError.A(1);
    end
    if abs(m(2))>data.MagDampingError.A(2)
        %Adjust to maximum
        m(2)=sign(m(2))*data.MagDampingError.A(2);
    end
    if abs(m(3))>data.MagDampingError.A(3)
        %Adjust to maximum
        m(3)=sign(m(3))*data.MagDampingError.A(3);
    end
    
    %Compute possible torque
    T_pos(end+1,1:3)=cross(m,mfield_round);
    %Compute difference to desired torque
    T_diff(end+1)=norm(T_pos(end,:)'-Tdes);
    %Save actution levels
    m_pos(end+1,1:3)= round(m(:)./data.MagDampingError.MTR(:)).*data.MagDampingError.MTR(:)+data.MagDampingError.MTB(:); 
    
end

%Get the configuration with minimum difference
[~,I]=min(T_diff);
T=cross(m_pos(I,:),mfield(:))';

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end