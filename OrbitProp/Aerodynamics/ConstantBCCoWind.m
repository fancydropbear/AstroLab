function [A,T,E] = ConstantBCCoWind (t,y,data)
% ConstantBCCoWind is a model to compute the drag of the spacecraft only
% using its ballistic coefficient BC = mass/(Cd*A).
% This model assumes that the atmosphere co-rotates with the Earth and 
% allows for the possibility to include wind.
%
% Parameters:
%   data.ConstantBCCoWind.BC -> Spacecarft ballistic coefficient BC = mass/(Cd*A)
%
%   data.ConstantBCCoWind.atmos_model -> Function handle specifying the
%   atmospheric model. The output of the atmospheric model must comply with
%   the standard format and at least contain the atmos.rho field specifying
%   the atmospheric density in kg/m3.
%
%   data.ConstantBCCoWind.wind_model -> Fucntion handle specifying the wind
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

%Flow is velocity - atmosphere rotatio - wind
wind = data.ConstantBCCoWind.wind_model(t,y,data); %Wind
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
corot = cross([0;0;we],y(1:3)); %Atmosphere rotation
flow = y(4:6) - corot - wind;

%Compute the atmospheric properties
atmos = data.ConstantBCCoWind.atmos_model(t,y,data);

%Compute drag acceleration
q = 1/2*atmos.rho*norm(flow)^2; %Dynamic pressure
A = -q*flow/norm(flow)/data.ConstantBCCoWind.BC;

%Torque (this model doesn't create any torque)
T=[0,0,0]';

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end