function [eclipse] = check_eclipse(sc,su)
% The function check_eclipse check if spacecarft is under eclipse.
% 
% The inputs are:
%   sc -> The spacecraft vector in ECEF coordinates.
%   su -> The sun vector in ECEF coordinates.
%
% And the ouputs:
%   eclipse -> True if the spacecarft is in eclipse.
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

%Coordinate transformation
x=su/norm(su);
z=cross(x,sc);
z=z/norm(z);
y=cross(z,x);

%Get the rotatation matrix
DCM = [x;y;z];

%Transform the sc vector to the new system
sc_p = DCM*sc';

%Earth parameters
Re = 6378136.49; %Equatorial Earth radius m [source: SMAD 3rd edition]

%Check if sc is in eclipse
if sc_p(2)<Re && sc_p(1)<0
    eclipse=true;
    return
else
    eclipse=false;
    return
end

end