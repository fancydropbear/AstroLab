function [A,T,E] = GravEGM2008(t,y,data)
% Implements the EGM2008 gravity model.
%
% Requires the gravitysphericalharmonic MATLAB function.

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

%Take into account the rotation of the earth
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
theta = t*we;
Rot = [cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0 0 1];

r = Rot*y(1:3);

%Compute the gravity acceleration.
if isfield(data,'gravity')
    if isfield(data.gravity,'degree')
        [a(1,1), a(2,1), a(3,1)] = gravitysphericalharmonic( r', 'EGM2008', data.gravity.degree);
    end
else
    [a(1,1), a(2,1), a(3,1)] = gravitysphericalharmonic( r', 'EGM2008');
end

%Get accelaration back to original coordinates
A = Rot\a;

%Torque (this model doesn't create any torque)
T=[0,0,0]';

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end