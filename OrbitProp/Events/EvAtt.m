function [value,isterminal,direction] = EvAtt(t,y,data)
% EvAtt is an event handler that stops integration when the spacecarft 
% achieves a certain attitude w.r.t the flow. This event is used to control
% the spin when DDsat counter-rotates its panels.
%
% Parameters:
%   data.EvAtt.wind_model -> Function handle specifying the wind
%   model. The output must be a 3x1 vector containg the wind in PCPF
%   coordinates.
%
%   data.EvAtt.r1-> Target roll angle.
%
%   data.EvAtt.r2 -> Target pitch angle.
%
%   data.EvAtt.r2 -> Target yaw angle.
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

%--- Attitude in flow reference frame ---%
%Flow is velocity - atmosphere rotatio - wind
wind = data.EvAtt.wind_model(t,y,data); %Wind
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

%--- Event values ---%
value = [r1-data.EvAtt.r1,r2-data.EvAtt.r2,r3-data.EvAtt.r3]; 
isterminal = [1,1,1]; %When the attitude is reached stop integration
direction = [0,0,0]; %Both directions


end