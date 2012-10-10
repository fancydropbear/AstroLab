function [t,x]=OrbitProp(x0,t0,tf,data)
% General purpose near-Earth modular orbit propagator.
%
% Usage:
%   [t,x] = OrbitProp(x0,t0,tf,data)
%
% Where:
%   Output:
%       t -> Array of time in s
%       x -> Array of state vectors (see below)
%   Input:
%       x0 -> Initial state vector (see below)
%       t0 -> Initial time in s
%       tf -> Final time in s (set it to inf for infinite)
%       data -> Structure array that configures the propagator and passes
%       on different paremeters. 
%
%
% State vectors have the following format:
%       x = [x,y,z,Vx,Vy,Vz,q1,q2,q3,q4,w1,w2,w3]
%
% Where:
%       xyz -> Position in m on a ECEF frame.
%       Vxyz -> Velocity in m/s.
%       q1234 -> Attitude quaterion. The reference frame is ECEF and the
%       Euler sequence is E123. q1 is the scalar value. 
%       w123 -> Angular rate in rad/s in the body reference frame
%
% Attitude can be ignored passing only the linear parameters as:
% x = [x,y,z,Vx,Vy,Vz]
%
% The integration can be stopped through the use of events. At the moment
% only one event can be set at a time. The event is specified through a
% function handler in data.event.event_function_handler.
% For convenience there are certain events already built in.
%
% Built in events:
%     @EvAlt -> This event handler allows to stop the integration when a
%     certain altitude is reached. This altitude is spceified in m through
%     the data.event.h.
%
%     @EvAtt -> This event handler allows to stop the integration when a
%     certain attitude with respect to the flow is achieved. The roll, 
%     pitch and yaw limits are spcified using data.EvAtt.r1,r2 and r3. Also
%     the wind model shall be specified using data.EvAtt.wind_model.
%
% The propagator relies on enviroment models set by the user. The models
% are passed in data.models as array of function handlers.
%
% For convinience this propagator comes with several models already
% built-in. But users can also create their own models.
%
% Built in gravity models:
%     @GravSpherical -> Earth spherical gravity field
%     @GravJ4 -> Earth gravity with the harmonics up to J4
%     @GravEGM96 -> EGM96 gravity model. By default up to 70 degree. The
%                     degree can ve changed via data.GravityEGM96.degree.
%     @GravityEGM2008 -> EGM2008 gravity model. By default up to 120 degree. The
%                     degree can ve changed via data.GravityEGM2008.degree.
%
% Built in aerodynamic models:
%     @ConstantBCCoWind -> Models drag with a constant ballistic coeficient
%                    BC=m/(A*Cd). The BC in kg/m2 is passed as 
%                    data.ConstantBCCoWind.BC. This function needs an atmospheric
%                    model and a wind model to extract the density and wind
%                    profile. This models need to be passed as a function 
%                    handlers in data.ConstantBCCoWind.atmos_model and
%                    data.ConstantBCCoWind.wind_model (see below for 
%                    built-in atmospheric models).
%
%     @AeroFDBCoWind -> Forces are computed looking up force coefficients in a 
%                      database. Coeficcients are assumed to be constant 
%                      through all altitudes. The coeficients are saved in 
%                      matrices in mat files. Each file contains the
%                      coeficcients for each direction (x, y and 
%                      z). Each file is organized is a nx4 matrix where a
%                      row is [angle, x Cf*A, y Cf*A, 
%                      z Cf*A]. The databases must be included
%                      in data.AeroFDBCoWind.rolldb_path, 
%                      data.AeroFDBCoWind.pitchdb_path and
%                      data.AeroFDBCoWind.yawdb_path. Also the mass must be
%                      provided in data.AeroFDBCoWind.mass in kg and the
%                      Cf*A for an unperturbed state provided in
%                      data.AeroFDBCoWind.CfA0. As usual atmospheric and
%                      wind models shall be provided in
%                      data.AeroFDBCoWind.atmos_model and
%                      data.AeroFDBCoWind.wind_model.
%
%     @AeroTDBCoWind -> Torques are computed looking up torque coefficients in a 
%                      database. Coeficcients are assumed to be constant 
%                      through all altitudes. The coeficients are saved in 
%                      matrices in mat files. Each file contains the
%                      coeficcients for each direction (x, y and 
%                      z). Each file is organized is a nx4 matrix where a
%                      row is [angle, roll Cf*A*l, pitch Cf*A*l, 
%                      yaw Cf*A*l]. The databases must be included
%                      in data.AeroTDBCoWind.rolldb_path, 
%                      data.AeroTDBCoWind.pitchdb_path and
%                      data.AeroTDBCoWind.yawdb_path. Also the inertia must be
%                      provided in data.sc_prop.I in kg/m2 and the
%                      Cf*A*l for an unperturbed state provided in
%                      data.AeroTDBCoWind.CfA0l0. As usual atmospheric and
%                      wind models shall be provided in
%                      data.AeroFDBCoWind.atmos_model and
%                      data.AeroFDBCoWind.wind_model.
%
% Built in atmospheric models:
%     @nrlmsise00 -> NRMSISE-00 atmospheric model.
%                    Requires: data.nrlmsise00.day0 as the intial day of the
%                    year, data.nrlmsise00.solar_activity to select different
%                    solar activity scenarios (options: 'Low', 'Medium',
%                    'High', 'HighShort' or 'Custom'. If 'Custom' is
%                    selected then data.nrlmsise00.f107Average,
%                    data.nrlmsise00.f107Daily and data.magneticIndex must be
%                    specified). More details can be found in the extended
%                    documentation.
%
%     @DbAtmosT -> The atmospheric properties are selected from a database.
%                  The databse contained data.DbAtmosT.Db is mx2 matrix. 
%                  The first column contains the time in seconds and the 
%                  second the density in kg/m3.
%
%     @WindHWM93 -> WindHWM93 implements the HWM93 wind models and gives 
%                   the output in a 3x1 wind vector in PCPF coordinates. It
%                   needs data.WindHWM93.day0 as the intial day of the year
%                   and data.WindHWM93.solar_activity to select different 
%                   solar activity scenarios (options: 'Low', 'Medium', 
%                   'High', 'HighShort' or 'Custom'. If 'Custom' is 
%                   selected then data.WindHWM93.f107Average,
%                   data.WindHWM93.f107Daily, data.WindHWM93.ap_daily and the
%                   data.WindHWM93.ap_3hour must be specified).
%                   WARNING! This model needs compilation of Fortran code.
%                   Please make sure that the compiled version for your
%                   platform exists. If not refer to the source code and
%                   compile.
%     
%     @NoWind -> NoWind models ignores wind and hence should be used when 
%                wind wants to be ignored but the models used require wind 
%                information.
%
% Built in solar radiation pressure models:
%      @SrpAccel -> SrpAccel is a model to compute the accelerations 
%                   produced by the solar radiation pressure (SRP). This 
%                   forces are dependant on the spacecarft attitude w.r.t
%                   the sun. To detrmine this forces the SRP ballistic 
%                   coefficients (m/r/A) are passed in a 3x3 matrix. It
%                   requires data.SrpAccel.sun_vector as a vector pointing 
%                   to the sun in ECEF coordinates, and also, 
%                   data.SrpAccel.invSrpBC as 3x3 matrix containing the 
%                   inverse SRP ballistic coefficient r*A/m. Each element 
%                   of the matrix x,y represents the force in x when the 
%                   sun is illuminating from y.
%
%     @SrpTorque -> SrpTorque is a model to compute the torques produced by
%                   the solar radiation pressure (SRP). This torques are 
%                   dependant on the spacecarft attitude w.r.t the sun. 
%                   To detrmine this torques the A*l*r coefficients are 
%                   passed in a 3x3 matrix. It requires 
%                   data.SrpAccel.sun_vector as a vector pointing to the 
%                   sun in ECEF coordinates, and also, data.SrpTorque.Alr 
%                   as a 3x3 matrix containing the inverse the A*L*r 
%                   product. Each element of the matrix x,y represents the 
%                   torque in x when the sun is illuminating from y.
%
% There are other models built in. These ones are located under the Others
% folder. Check their documentation to know what they are for and learn how
% to use them.
%
% Depedencies of built in models:
%     @GravityJ4 -> gravityzonal (Aerospace Toolbox)
%     @GravityEGM96 -> gravitysphericalharmonic (Aerospace Toolbox)
%     @GravityEGM2008 -> gravitysphericalharmonic (Aerospace Toolbox)
%     @nrlmsise00 -> atmosnrlmsise00 (Aerospace Toolbox)
%
% To check if your installation has this function you can use 'exist(fun)'
%
% Examples are provided in the tests folder.

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

%---- Check initial state vector ---%
%Convert intial state to column vector
x0=x0(:);
%Check size of initial state column vector
if size(x0,1)<6
    %Initial state vector doesn't contain minimum number of elements
    %Throw error
    error('x0 must be have at least x,y,z,Vx,Vy,Vz');
elseif size(x0,1)==6
    %Ignoring attitude
    x0(7:10)=[1,0,0,0]'; %Arbitrary initial attitude
    x0(11:13)=[0,0,0]'; %Zero attitude rates
    data.sc_prop.I=diag([1,1,1]); %Arbitrary inertia
elseif size(x0,1)<13
    %The size of the initial state vector is incorrect.
    %Throw error
    error('The initial vector size doesn''t match with x,y,z,Vx,Vy,Vz or x,y,z,Vx,Vy,Vz,q1,q2,q3,q4,w1,w2,w3'); 
end

%--- Check for models ---%
try isfield(data,'models');
    %The field exists
    %Conver to column cell array
    data.models=data.models(:);
catch
    %There are no models specified
    %Throw error
    error('No models specified in data.models');
end

%--- Check for events ---%
try isfield(data.event,'function');
    %Event function found
    options=odeset('RelTol',1e-6,'AbsTol',1e-6,'Event',data.event.function);
catch
    %Event function not found
    options=odeset('RelTol',1e-6,'AbsTol',1e-6);
end

%--- Integration ---%
[t,x]=ode45(@DyEq,[t0,tf],x0,options,data);

end