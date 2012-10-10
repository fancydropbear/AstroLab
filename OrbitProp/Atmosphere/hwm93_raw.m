% The HWM93 function calls the Horizontal Wind Model 1993 FORTRAN code via 
% mex gateway
%
% Syntax:
%
%   W = HWM07(IYD,SEC,ALT,GLAT,GLON,SLT,F107A,F107,AP_DAILY,AP_3HOUR)
%
% Inputs:
%
%   IYD - date in yyddd format
%   SEC - UT [sec]
%   ALT - altitude [km]
%   GLAT - North Geographic Latitude [degrees]
%   GLON - East Geographic Longitude [degrees]
%   F107A    - 107a 81-dy average solar flux (W/cm^2/Hz)
%   F107     - f107 solar flux (W/cm^2/Hz)
%   AP_DAILY - daily Ap index
%   AP_3HOUR - 3 hr Ap index
%
% Outputs:
%   
%   W - 1x2 vector containing meridional [+north] (1) and zonal [+east] (2) wind [m/s]
%
%
% 
% HWM93 mex gateway written by:
%
% Tim Duly (duly2@illinois.edu)
% Remote Sensing and Space Sciences Group (http://rsss.csl.uiuc.edu/)
% Department of Electrical and Computer Engineering
% University of Illinois at Urbana-Champaign
% Current:   2010 Aug 6 / Version 1.0
% Original:  2010 Aug 6 / Version 1.0
