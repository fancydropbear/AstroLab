function [A,T,E] = GravSpherical(t,y,data)
% GravSpherical implements a spherical gravity fiels with GM = 398600.441e9
% m3/s2
%
% This model don't need any configuration.

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

%Compute the ditsance from origin
r=norm(y(1:3));

%Gravity constant
mu=398600.441e9;  %GM Earth m3/s2 [source: SMAD 3rd edition]

%compute local gravity
g=mu/r^2;

%Fromat result
A = -g*y(1:3)/r;

%Torque (this model doesn't create any torque)
T=[0,0,0]';

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end