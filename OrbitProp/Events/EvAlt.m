function [value,isterminal,direction] = EvAlt(t,y,data)
% EvAlt is an event handler that stops the integration when a certain
% altitude is reached.
%
% Parameters:
%   data.event.h -> Altitude in m to stop integration.

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

% Earth constants
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]

% Event values
value = norm(y(1:3)) - data.event.h - Req ; 
isterminal = 1; %When the altitude is reached stop integration
direction = 0; %Both directions



end