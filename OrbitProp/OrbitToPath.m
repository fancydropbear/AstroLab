% Script that adds all necesarry functions to the path for OrbitProp to
% work.
%
%
% This script is robust. It doesn't depends on the your current folder.
% This script should be executed before executing OrbitProp. If it not
% executed OrbitProp will fail displaying an error related to not finding
% the necessary functions.
%
% Warning! Errors will occurr if in the path to this file there are folders
% with the name 'OrbitToPath'.

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

%--- CODE ---%

%Get current full path
currentpath=strrep(mfilename('fullpath'),mfilename(),'');

%Add all the other folders to path
addpath(fullfile(currentpath,'Gravity')); %Adds the gravity models
addpath(fullfile(currentpath,'Atmosphere')); %Adds the atmospheric models
addpath(fullfile(currentpath,'Aerodynamics')); %Adds the aerodynamics models
addpath(fullfile(currentpath,'Dynamics')); %Adds dynamic equations
addpath(fullfile(currentpath,'SRP')); %Adds the solar radiation pressure models
addpath(fullfile(currentpath,'Others')); %Adds other miscellanious models
addpath(fullfile(currentpath,'Events')); %Adds the event handlers
addpath(fullfile(currentpath,'Custom')); %Adds the custom  models


