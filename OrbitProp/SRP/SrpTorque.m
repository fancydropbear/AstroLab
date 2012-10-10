function [A,T,E] = SrpTorque(t,y,data)
% SrpTorque is a model to compute the torques produced by the solar
% radiation pressure (SRP). This torques are dependant on the spacecarft
% attitude w.r.t the sun. To detrmine this torques the A*l*r coefficients
% are passed in a 3x3 matrix.
%
% Parameters:
%   data.SrpAccel.sun_vector -> Vector pointing to the sun in ECEF
%   coordinates.
%
%   data.SrpTorque.Alr -> 3x3 matrix containing the inverse the A*L*r 
%   product. Each element of the matrix x,y represents the Torque in x when
%   the sun is illuminating from y.
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

%Check if spacecraft is in eclipse
if check_eclipse(y(1:3),data.SrpTorque.sun_vector)
    %Spacecraft in eclipse
    T = [0;0;0];
else
    %Spacecraft is not in eclipse
    Sun_P = 9.15/1e6; %Sun pressure at 1 AU in N/m2
    %Compute torque
    T=data.SrpTorque.Alr*Sun_P;
    
    %Get sun vector in body axis
    DCM = quat2dcm([y(10),y(7:9)']);
    sun = DCM*data.SrpTorque.sun_vector;
    
    %Get the correct sign and illumation fractions
    T=T*sun/norm(sun);
    
end

%Acceleration (this model don't produce linear acceleration)
A=[0;0;0];

%Extra state variables (this model doesn't need extra state variables)
E=zeros(1,length(y)-13);

end