function [A,T,E] = GravityGradient(t,y,data)
% GravityGradient is a model to compute the torques produced by the gravity
% gradient. This torques are dependant on the spacecarft attitude with
% respect to the Earth.

%--- Copyright notice ---%
% Copyright 2013 Cranfield University
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

%Parameters
mu=398600.441e9;  %GM Earth m3/s2 [source: SMAD 3rd edition]

% Get position vector in body axes
DCM = quat2dcm([y(10),y(7:9)']);
R = DCM*y(1:3);

%Compute torque in body axes
T=3*mu/norm(R)^3*cross(R/norm(R),data.sc_prop.I*(R/norm(R)));

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);
end