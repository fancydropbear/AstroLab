function [A,T,E]=magHysteresis(t,y,data)
% magHysteresis is a model that computes the magnetic hysteresis of
% magnetic permeable material.
% To model the hysteresis 3 new extra state variables have to be intoduced.

% In this case are the flux density induced in the rods. This model then
% computes the torque produced by the hysteresis loops, but also the
% derviates to these extra state variables. The location of these extra
% state variables in the state vector y is specified in
% data.magHysteresis.BIndex. This variables must be initialized.
%
% The Earth magnetic field is modeled using the IGRF-11 and therfore the
% model requires the igrf11magm MATLAB function.
%
% Parameters:
%   data.MagDamping.Y -> Year (include decimal point for increased precision)
%   
%   data.Hysteresis.Bs -> A 3x1 vector with the rods saturation level in 
%   Teslas.
%
%   data.Hysteresis.Br -> A 3x1 vector with the rods remanences in Teslas.
%   
%   data.Hysteresis.Hc -> A 3x1 vector with the rods coercivity in A/m.
%
%   data.Hysteresis.V -> A 3x1 vector with the volume in m3 of the torquer
%   rods in all three directions.
%
%   data.magHysteresis.BIndex -> A vector containing the 3 positions of the
%   extra state variables (that is the flux density induced in the rods).
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

%Earth constants
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
%Vacuum permeability
mu0=4*pi*1e-7;
%Rot from ECEF to body
Rbe=quat2dcm(y(7:10)');

%Take into account the rotation of the earth
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
theta = t*we;
Rot = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];


%--- Magnetic field at t ---%
lla=ecef2lla(y(1:3)',f,Req); %Compute taking assuming Earth as an ellipsoid
lat1=lla(1);
lon1=lla(2);
h1=lla(3);
%Take into account rotation of the earth
lon1 = mod(lon1 - t*we*180/pi,360);
if lon1>180
    lon1 = lon1-360;
elseif lon1<-180
    lon1 = lon1+360;
end
%Magnetic field in T
H_aux =igrf11magm(h1, lat1, lon1, data.magHysteresis.Y)'/1e9/mu0;
H1 = Rbe*Rot\[-H_aux(1)*sind(lat1)-H_aux(3)*cosd(lat1);H_aux(2);H_aux(1)*cosd(lat1)-H_aux(3)*sind(lat1)];

%--- Magnetic field at t+1s ---%
lla=ecef2lla(y(1:3)'+y(4:6)',f,Req); %Future position in 1 s
lat2=lla(1);
lon2=lla(2);
h2=lla(3);
%Take into account rotation of the earth
lon2 = mod(lon2 - t*we*180/pi,360);
if lon2>180
    lon2 = lon2-360;
elseif lon2<-180
    lon2 = lon2+360;
end
%Magnetic field in T
H_aux =igrf11magm(h2, lat2, lon2, data.magHysteresis.Y)'/1e9/mu0;
H2 = Rbe*Rot\[-H_aux(1)*sind(lat2)-H_aux(3)*cosd(lat2);H_aux(2);H_aux(1)*cosd(lat2)-H_aux(3)*sind(lat2)];

%Hdot
Hdot=H2-H1-cross(y(11:13),H1);

%--- Hysteresis equations ---%
B=y(data.magHysteresis.BIndex+13);
k=1./data.magHysteresis.Hc.*tan(pi/2.*data.magHysteresis.Br./data.magHysteresis.Bs); %Shaping factor
if Hdot>0
    Hlim=1./k.*tan(pi./2.*B./data.magHysteresis.Bs)-data.magHysteresis.Hc;
else
    Hlim=1./k.*tan(pi./2.*B./data.magHysteresis.Bs)+data.magHysteresis.Hc;
end

Bdot = 2/pi*k.*data.magHysteresis.Bs.*cos(pi/2.*B./data.magHysteresis.Bs).^2.*((H1(:)-Hlim)./(2*data.magHysteresis.Hc)).^2.*Hdot(:);

%--- Torque ---%
%Rod magnetic dipole
D=B.*data.magHysteresis.V;
%Torque
T=cross(D,H1);

%--- B derivate ---%
E=zeros(1,length(y)-13);
E(data.magHysteresis.BIndex)=Bdot; %Passed as Extra state variable

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

end