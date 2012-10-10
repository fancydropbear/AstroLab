function [atmos] = DbAtmosT (t,y,data)
% DbAtmosT is an atmospheric model where properties are extracted from a
% DB. This model is used when actual data about the atmosphere is known.
%
% Parameters:
%   data.DbAtmosT.Db -> mx2 matrix. The first column contains the time in 
%   seconds and the second the density in kg/m3.
%

%--- Copyright notice ---%
% Copyright 2012 Cranfield University
% Written by Josep Virgili
%
% This file is part of the AstroLab toolkit
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

%---- CODE ----%

% Interpolate inside the matrix to extract density.
atmos.rho = interp1(data.DbAtmosT.Db(:,1),data.DbAtmosT.Db(:,2),t);

end