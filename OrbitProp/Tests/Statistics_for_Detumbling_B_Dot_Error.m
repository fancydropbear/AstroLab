% Statics for Bdot Detumbling with all error sources is a test script for OrbitProp

% Statics_for_the_Detumbling_B_Dot_Error is a test script for OrbitProp
% This function is dedicated to get the (detumbling time statistics), (power
% statics)
% Before running this static file, remove the 'clear function' from the
% DDsatDetumblingError.m

%--- Copyright notice ---%
% Copyright 2012-2013 Cranfield University
% Written by Josep Virgili and Daniel Zhou Hao
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
clear all;
clc;

% 30 times is approved by Chi-Squared confidence interval
for counter = 1:30
    
    % Find where the time for angular rate getting to small enough
    DDsatDetumblingBDotError;
    
    for counter2 = 1: length(t)
        if norm(rad2deg(x(counter2,11:13))) < 0.5
            detumb_t(counter) = t(counter2);
            avg_power(counter) = mean(PT(1:counter2));
            break;
        end
    end
end

mean_t = mean(detumb_t(:))/60
std_t = std(detumb_t(:))/60
mean_avg_power = mean(avg_power(:))
std_avg_power = std(avg_power(:))
