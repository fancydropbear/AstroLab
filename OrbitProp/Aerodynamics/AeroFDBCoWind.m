function [A,T,E] = AeroFDBCoWind(t,y,data)
% AeroFDBCoWind is a model to compute the aerodynamic forces that act on
% the spacecraft. This forces can be dependant on the spacecarft attitude
% with respect to the flow. To detrmine this relation the Cf*A are passed
% as aerodynamic databases.
% This model assumes that the atmosphere co-rotates with the Earth and 
% allows for the possibility to include wind.
%
% Parameters:
%   data.AeroFDBCoWind.rolldb -> mx4 matrix. The first column contains
%   the roll angle in rad. The 2-4 column contain the x,y,z Cf*A that will 
%   act in that direction. The values must be absolute.
%
%   data.AeroFDBCoWind.pitchdb -> mx4 matrix. The first column contains
%   the pitch angle in rad. The 2-4 column contain the x,y,z Cf*A that will 
%   act in that direction. The values must be absolute.
%
%   data.AeroFDBCoWind.yawdb -> mx4 matrix. The first column contains
%   the yaw angle in rad. The 2-4 column contain the x,y,z Cf*A that will 
%   act in that direction. The values must be absolute.
%
%   data.AeroFDBCoWind.mass -> Mass in kg of the spacecarft.
%
%   data.AeroFDBCoWind.CfA0 -> 3x1 vector containing the Cf*A that acts in
%   x,y,z when the roll, picth and yaw are 0.
% 
%   data.AeroFDBCoWind.atmos_model -> Function handle specifying the
%   atmospheric model. The output of the atmospheric model must comply with
%   the standard format and at least contain the atmos.rho field specifying
%   the atmospheric density in kg/m3.
%
%   data.AeroFDBCoWind.wind_model -> Function handle specifying the wind
%   model. The output must be a 3x1 vector containg the wind in PCPF
%   coordinates.

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

%--- Attitude in flow reference frame ---%
%Flow is velocity - atmosphere rotatio - wind
wind = data.AeroFDBCoWind.wind_model(t,y,data); %Wind
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

%--- Database lookup ---%
rolldb=data.AeroFDBCoWind.rolldb;
pitchdb=data.AeroFDBCoWind.pitchdb;
yawdb=data.AeroFDBCoWind.yawdb;

%--- CdA ---%
CfA_roll = [interp1(rolldb(:,1),rolldb(:,2),r1);
             interp1(rolldb(:,1),rolldb(:,3),r1);
             interp1(rolldb(:,1),rolldb(:,4),r1)];
         
CfA_pitch =[interp1(pitchdb(:,1),pitchdb(:,2),r2);
             interp1(pitchdb(:,1),pitchdb(:,3),r2);
             interp1(pitchdb(:,1),pitchdb(:,4),r2)];
         
CfA_yaw =  [interp1(yawdb(:,1),yawdb(:,2),r3);
             interp1(yawdb(:,1),yawdb(:,3),r3);
             interp1(yawdb(:,1),yawdb(:,4),r3)];
         
%--- Compute Force ---%
%Compute the atmospheric properties
atmos = data.AeroFDBCoWind.atmos_model(t,y,data);

%Dynamic pressure
q=1/2*atmos.rho*norm(flow)^2;

%Acceleration
A=q/data.AeroFDBCoWind.mass*(CfA_roll+CfA_pitch+CfA_yaw-2*data.AeroFDBCoWind.CfA0);
A=Rfe'*A; %Change coordinates to PCPF

%Torque (this model doesn't create any torque)
T=[0,0,0]';

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end