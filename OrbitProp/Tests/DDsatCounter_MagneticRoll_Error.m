% DDsatCounter_MagneticRoll_Error is a test script for OrbitProp
%
% In this case this example shows the orbit propagator in it's full
% potential using several models and simulating DDsat with the panels
% counter-rotating. 
%
% In this case the roll is controlled through the magnetic torquers. And
% all error/bias/uncertainties sources has beenn added into this sinulation
%
% The aerodynamic databases are in the AeroDB folder
%

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

%--- CODE ---%

%Clen up
clc
clear all
close all

%Add paths
addpath('..')
OrbitToPath

%--- Initial state vector and model configuration ---%
%Earth parameters
Re = 6378136.49; %Equatorial Earth radius m [source: SMAD 3rd edition]
mu=398600.441e9;  %GM Earth m3/s2 [source: SMAD 3rd edition]
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]

%Initial conditions
h = 320*1e3; %Initial altitude in m
i = 79; %Inclination
v0 = sqrt(mu/(Re+h)); %Initial velocity
x0 = [Re+h,0,0,0,v0*cosd(i),v0*sind(i)]; %Initial vector state
%Atmosphere co-rotation velocity
Vco = norm(cross([0;0;we],x0(1:3))); 

tf =[0,2*90*60]; %Integration time
rx=deg2rad(-90);
ry=deg2rad(-90);
rz=deg2rad(-atand((v0*sind(i))/(v0*cosd(i)-Vco)));
q0=angle2quat(rx,ry,rz,'XYZ'); %Initial attitude quaterion
w0=[0,deg2rad(-360/90/60),0]; %Initial angular rate
%Format the initial 
x0=[x0,q0,w0];

%Spacecarft propetries
% data.sc_prop.I=[0.0033,0,0;
%                 0,0.0083,0;
%                 0,0,0.0083]; %Old Inertia

data.sc_prop.I=[0.04,0,0;
                0,0.0177,0;
                0,0,0.0177]; %Updated Inertia

            
%Aerodynamic databases
rolldb_FP = 'AeroDB/Counter/pmz-45pmy0/F_roll_pmz-45_pmy0.mat';
pitchdb_FP = 'AeroDB/Counter/pmz-45pmy0/F_pitch_pmz-45_pmy0.mat';
yawdb_FP = 'AeroDB/Counter/pmz-45pmy0/F_yaw_pmz-45_pmy0.mat';
rolldb_FM = 'AeroDB/Counter/pmz45pmy0/F_roll_pmz45_pmy0.mat';
pitchdb_FM = 'AeroDB/Counter/pmz45pmy0/F_pitch_pmz45_pmy0.mat';
yawdb_FM = 'AeroDB/Counter/pmz45pmy0/F_yaw_pmz45_pmy0.mat';
rolldb_TP = 'AeroDB/Counter/pmz-45pmy0/T_roll_pmz-45_pmy0.mat';
pitchdb_TP = 'AeroDB/Counter/pmz-45pmy0/T_pitch_pmz-45_pmy0.mat';
yawdb_TP =  'AeroDB/Counter/pmz-45pmy0/T_yaw_pmz-45_pmy0.mat';
rolldb_TM = 'AeroDB/Counter/pmz45pmy0/T_roll_pmz45_pmy0.mat';
pitchdb_TM = 'AeroDB/Counter/pmz45pmy0/T_pitch_pmz45_pmy0.mat';
yawdb_TM = 'AeroDB/Counter/pmz45pmy0/T_yaw_pmz45_pmy0.mat';
            
%Models
data.models={@AeroFDBCoWind};
data.models(end+1)={@GravJ4};
data.models(end+1)={@SrpAccel};
data.models(end+1)={@AeroTDBCoWind};
data.models(end+1)={@MagDampingError};
data.models(end+1)={@SrpTorque};
data.models(end+1)={@MagRollControl};
data.models(end+1)={@GravityGradient};

%Configure AeroFDBCoWind
FP_rolldb=load(rolldb_FP);
FP_pitchdb=load(pitchdb_FP);
FP_yawdb=load(yawdb_FP);
data.AeroFDBCoWind.rolldb=FP_rolldb.F_roll;
data.AeroFDBCoWind.pitchdb=FP_pitchdb.F_pitch;
data.AeroFDBCoWind.yawdb=FP_yawdb.F_yaw;
data.AeroFDBCoWind.atmos_model=@nrlmsise00;
data.AeroFDBCoWind.wind_model=@NoWind;
data.AeroFDBCoWind.mass=2;
data.AeroFDBCoWind.CfA0=[FP_rolldb.F_roll(1,2);FP_rolldb.F_roll(1,3);FP_rolldb.F_roll(1,4)];

%Configure nrlmsise00
data.nrlmsise00.solar_activity = 'High'; %Level of solar activity
data.nrlmsise00.day0 = 10; %Initial day of the year

%Configure SrpAccel
data.SrpAccel.sun_vector = [1;0;0];
SrpBCx= 2/(0.01*0.9);
SrpBCyz= 2/(0.1*0.9);
data.SrpAccel.invSrpBC = diag([1/SrpBCx,1/SrpBCyz,1/SrpBCyz]);

%Configure AeroTDBCoWind model
TP_rolldb=load(rolldb_TP);
TP_pitchdb=load(pitchdb_TP);
TP_yawdb=load(yawdb_TP);
data.AeroTDBCoWind.rolldb = TP_rolldb.T_roll;
data.AeroTDBCoWind.pitchdb = TP_pitchdb.T_pitch;
data.AeroTDBCoWind.yawdb = TP_yawdb.T_yaw;
data.AeroTDBCoWind.atmos_model=@nrlmsise00;
data.AeroTDBCoWind.wind_model=@NoWind;
data.AeroTDBCoWind.CfA0l0=[TP_rolldb.T_roll(1,2);TP_rolldb.T_roll(1,3);TP_rolldb.T_roll(1,4)];

%Configure WindHWM93
data.WindHWM93.day0 = 10; %Initial day of the year
data.WindHWM93.solar_activity = 'Medium'; %Level of solar activity

%Configure MagDampingError
data.MagDampingError.c = 2e-4;
data.MagDampingError.A = [0.2,0.2,0.2];
data.MagDampingError.Y = 2012;
data.MagDampingError.w = [0,deg2rad(-360/90/60),0];
data.MagDampingError.verb = 0; 
%Configure Error Sources
% data.MagDampingError.MTB -> Actuation Level Bias Vector in Am2
% data.MagDampingError.MMR -> Magnetometre Resolution in Teslas 
% data.MagDampingError.MTR -> Magnetorquer Resolution in Am2 
% data.MagDampingError.GE -> Gyroscope Reading One Sigma Uncertaincy
% Vector in rad/s

data.MagDampingError.GE = deg2rad ([0.1,0.1,0.1]');   
% 1st Case[0.1,0.1,0.1]'
% 2nd Case[0.05,0.05,0.05]'
% 3rd Case[0.025,0.025,0.025']
%Check the sensitivity with Gyroscope uncertainty reading
%Start with 0.1 reg sigma.
data.MagDampingError.MMR =  ones(3,1)*0.7e-6;               % 7 mili gauss = 0.7 microteslas
data.MagDampingError.MTR = ones(3,1)*0.2/(2^8);   % +- 8 bits of 2 Am2
data.MagDampingError.MTB = ones(3,1)*normrnd (0, 0.025*0.2);
data.MagDampingError.SF = 0.1; % 10 Hz


%Configure SrpTorque
data.SrpTorque.sun_vector = [1;0;0];
Alr= 0.1*0.9*0.1-0.01*0.9*0.1;
data.SrpTorque.Alr = [0,0,0;
                     0,0,Alr;
                     0,Alr,0];
                 
%Configure MagRollControlError
data.MagRollControl.A = [0.2,0.2,0.2];
data.MagRollControl.Y = 2012;
data.MagRollControl.k = 2e-4;
data.MagRollControl.wind_model=@NoWind;
data.MagRollControl.verb = 0;

%--- Integrate ---%
[t,x] = OrbitProp(x0,tf(1),tf(2:end),data);

%--- POST-PROCESS ---%

%--- Ground Track plot ---%
%Compute altitude lat and long
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
lla(:,2) = lla(:,2)- t*we*180/pi;
%plot ground track
figure
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
axis tight
geoshow(lla(:,1),lla(:,2));
title('Ground track')

%--- Euler Angles plot ---%
%Get euler angles
[r1 r2 r3] = quat2angle([x(:,7:10)],'XYZ');
%Plot euler angles
figure
plot(t/60,rad2deg(r1),t/60,rad2deg(r2),t/60,rad2deg(r3));
xlabel('Time [min]')
ylabel('Euler angles [deg]')
legend('E1','E2','E3')
title('Euler angles')

%--- Rate plot ---%
figure
plot(t/60,rad2deg(x(:,11)),t/60,rad2deg(x(:,12)),t/60,rad2deg(x(:,13)));
xlabel('Time [min]')
ylabel('Angular rates [deg/s]')
legend('Roll','Pitch','Yaw')
title('Angular rates')

%-- Flow attitude plot --%
for i=1:size(x,1)
    %Flow reference frame
    wind = data.AeroTDBCoWind.wind_model(t(i),x(i,:)',data); %Wind
    corot = cross([0;0;we],x(i,1:3)); %Atmosphere rotation
    flow(i,1:3) = x(i,4:6)-corot-wind';
    xa = flow(i,1:3)/norm(flow(i,1:3)); %Velocity (roll axis)
    zaux = -x(i,1:3)/norm(x(i,1:3)); %Nadir pointing
    ya = cross(zaux,xa); % Velocity anti-paralell (pitch axis)
    za = cross(xa,ya); % Near nadir (yaw axis)

    %Rotation matrix from ECEF to flow ref
    Rfe=[xa;ya;za];

    %Rot from ECEF to body
    Rbe=quat2dcm([x(i,7:10)]);

    %Rot from flow to body
    Rbf = Rbe*inv(Rfe);
    %Get angles
    [r1(i), r2(i), r3(i)] = dcm2angle(Rbf,'XYZ');
    
    %Check for eclipse
    eclipse(i) = check_eclipse(x(i,1:3)',data.SrpTorque.sun_vector);
end

%--- Plot flow attitude ---%
figure
plot(t/60,rad2deg(r1),t/60,rad2deg(r2),t/60,rad2deg(r3),t/60,eclipse,'k');
xlabel('Time [min]')
ylabel('Euler angles [deg]')
legend('Roll','Pitch','Yaw','Eclipse')
title('Attitude w.r.t the flow')


%Plot density
for i=1:size(x,1)
    %Compute the atmospheric properties
    atmos = data.AeroTDBCoWind.atmos_model(t(i),x(i,:)',data);
    rho(i)=atmos.rho;
    %Dynamic pressure
    q(i)=rho(i)*norm(flow(i,1:3))^2;
    
    %Flow without wind
    corot = cross([0;0;we],x(i,1:3)); %Atmosphere rotation
    flow_nowind(i,1:3) = x(i,4:6)-corot;
    
    %Observed rho
    rho_obs(i)=q(i)/norm(flow_nowind(i,1:3))^2;
     
end

%--- Rho vs observed rho comparison ---%
figure
plot(t/60,rho_obs,t/60,rho)
xlabel('Time [min]')
ylabel('Density [kg/m3]')
legend('Observed','Real')
title('Observed vs real density')


%--- Density variability ---%
figure
plot(t/60,rho_obs/rho_obs(1))
xlabel('Time [min]')
ylabel('Relative density')
title('Observed density variability')

%--- Rho vs observed rho comparison ---%
figure
plot(t/60,rho_obs./rho)
xlabel('Time [min]')
ylabel('Density ratio')
title('Observed / real density')

%--- Magnetic Torque ---%
%Calculate Torque
Torque=[];
for i=1:length(t)
 [A,Torque(end+1,1:3),E] = MagDampingError(t(i),x(i,:)',data);
end

%Torque plot
figure
plot(t/60,Torque(:,1),t/60,Torque(:,2),t/60,Torque(:,3))
xlabel('Time [mins]');
ylabel('Torque [Nm]');
legend('x Torquer','y Torquer','z Torquer');
title('Torque Plot');

%--- Power & Energy Plots ---%
Al = []; %Actuation level
PT = []; %Total Power
PP = []; %Partial Power
% Power per actuation level Am2 of actuators (parameter that needs
% adjusting to every case)
px = 0.57/0.2;
py = 0.2/0.2;
pz = 0.2/0.2;
p = [px,py,pz];

for i=1:length(t)
    %--- Get magnetic field ---%
    %Latitude and longitude
    we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
    Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
    f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
    lla=ecef2lla(x(i,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
    lat=lla(1);
    lon=lla(2);
    h=lla(3);
    %Take into account rotation of the earth
    lon = mod(lon - t(i)*we*180/pi,365);
    if lon>180
        lon = lon-360;
    elseif lon<-180
        lon = lon+360;
    end
    
    %Magnetic field in T
    mfield_ECEF = igrf11magm(h, lat, lon, data.MagDampingError.Y)/1e9;
    %Change magnetic field to body axes.
    DCM = quat2dcm([x(i,10),x(i,7:9)]);
    mfield = DCM*mfield_ECEF';
    
    %Preliminary actuation levels m_p
    m_p=(cross(mfield,Torque(i,:))/norm(mfield)^2)';
    %Check different scalar values
    k_max=(data.MagDampingError.A'-m_p)./mfield; %Scalar value for maximum actuation level on one axis
    k_min=(-data.MagDampingError.A'-m_p)./mfield; %Scalar value for minus maximum actuation level on one axis
    k_0=-m_p./mfield; %Scalar value for 0 actuation level on one axis
    
    %Scan through the different scalar values
    err=1e-3/100;
    m_v=[];
    p_v=[];
    for k=[k_max',k_min',k_0']
        m=m_p+k*mfield;
        
        if abs(m(1))>data.MagDampingError.A(1)*(1+err)
            %Solution is not real
            continue
        end
        if abs(m(2))>data.MagDampingError.A(2)*(1+err)
            %Solution is not real
            continue
        end
        if abs(m(3))>data.MagDampingError.A(3)*(1+err)
            %Solution is not real
            continue
        end
        
       m_v(end+1,1:3) = m;
       p_v(end+1) = sum(abs(m_v(end,:)).*p);
    end
           
   
   % Select the solution has the minmum total power
   [PT(end+1),I] = min(p_v);
   PP(end+1,1:3)= abs(m_v(I,1:3)).*p;
   Al(end+1,1:3) = m_v(I,1:3);
   
end

%Power plot
figure
plot(t/60,PT,t/60,PP(:,1),t/60,PP(:,2),t/60,PP(:,3))
xlabel('Time [mins]');
ylabel('Total Power [W]');
legend('Total Power','x power','y power','z power');
title('Power Plot');

%Actuation level plot
figure
plot(t/60,Al(:,1),t/60,Al(:,2),t/60,Al(:,3))
xlabel('Time [mins]');
ylabel('Actuation Level [Am2]');
legend('x actuation level','y actuation','z actuation');
title('Actuation Level Plot');


