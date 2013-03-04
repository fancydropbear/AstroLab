% DDsatCo is a test script for OrbitProp
%
% In this case this example shows the orbit propagator in it's full
% potential using several models and simulating DDsat with the panels
% co-rotating.
% The aerodynamic databases are in the AeroDB folder

%--- CODE ---%

%Clean up
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

tf = 0:60:5 * 90 * 60; %Integration time
rx=deg2rad(-90);
ry=deg2rad(-90);
rz=deg2rad(-atand((v0*sind(i))/(v0*cosd(i)-Vco)));
q0=angle2quat(rx,ry,rz,'XYZ'); %Initial attitude quaterion
w0=[0,deg2rad(-360/90/60),0]; %Initial angular rate
%Format the initial 
x0=[x0,q0,w0];

%Spacecarft propetries
data.sc_prop.I=[0.0033,0,0;
                0,0.0083,0;
                0,0,0.0083]; %Inertia
            
%Models
data.models={@GravJ4};
% data.models(end+1)={@ConstantBCCoWind};
data.models(end+1)={@AeroFDBCoWind};
data.models(end+1)={@SrpAccel};
data.models(end+1)={@AeroTDBCoWind};
data.models(end+1)={@MagDamping};
data.models(end+1)={@SrpTorque};
data.models(end+1)={@MagRollControl};
data.models(end+1)={@GravityGradient};

% %Configure ConstantBC model
% data.ConstantBCCoWind.BC=2/(0.0916*2.2);
% data.ConstantBCCoWind.atmos_model=@nrlmsise00;
% data.ConstantBCCoWind.wind_model=@WindHWM93;

%Configure AeroFDBCoWind
rolldb=load('AeroDB/Co/pmz90pmy0/F_roll_pmz90_pmy0.mat');
pitchdb=load('AeroDB/Co/pmz90pmy0/F_pitch_pmz90_pmy0.mat');
yawdb=load('AeroDB/Co/pmz90pmy0/F_yaw_pmz90_pmy0.mat');
data.AeroFDBCoWind.rolldb=rolldb.F_roll;
data.AeroFDBCoWind.pitchdb=pitchdb.F_pitch;
data.AeroFDBCoWind.yawdb=yawdb.F_yaw;
data.AeroFDBCoWind.atmos_model=@nrlmsise00;
data.AeroFDBCoWind.wind_model=@WindHWM93;
data.AeroFDBCoWind.mass=2;
data.AeroFDBCoWind.CfA0=[rolldb.F_roll(1,2);rolldb.F_roll(1,3);rolldb.F_roll(1,4)];

%Configure nrlmsise00
data.nrlmsise00.solar_activity = 'High'; %Level of solar activity
data.nrlmsise00.day0 = 10; %Initial day of the year

%Configure SrpAccel
data.SrpAccel.sun_vector = [1;0;0];
SrpBCx= 2/(0.01*0.9);
SrpBCyz= 2/(0.1*0.9);
data.SrpAccel.invSrpBC = diag([1/SrpBCx,1/SrpBCyz,1/SrpBCyz]);

%Configure AeroTDBCoWind model
rolldb=load('AeroDB/Co/pmz90pmy0/T_roll_pmz90_pmy0.mat');
pitchdb=load('AeroDB/Co/pmz90pmy0/T_pitch_pmz90_pmy0.mat');
yawdb=load('AeroDB/Co/pmz90pmy0/T_yaw_pmz90_pmy0.mat');
data.AeroTDBCoWind.rolldb = rolldb.T_roll;
data.AeroTDBCoWind.pitchdb = pitchdb.T_pitch;
data.AeroTDBCoWind.yawdb = yawdb.T_yaw;
data.AeroTDBCoWind.atmos_model=@nrlmsise00;
data.AeroTDBCoWind.wind_model=@WindHWM93;
data.AeroTDBCoWind.CfA0l0=[rolldb.T_roll(1,2);rolldb.T_roll(1,3);rolldb.T_roll(1,4)];

%Configure WindHWM93
data.WindHWM93.day0 = 10; %Initial day of the year
data.WindHWM93.solar_activity = 'Medium'; %Level of solar activity

%Configure MagDamping
data.MagDamping.c = 2e-4;
data.MagDamping.A = [0.2,0.2,0.2];
data.MagDamping.Y = 2012;
data.MagDamping.w = [0,deg2rad(-360/90/60),0];
data.MagDamping.verb=1;

%Configure SrpTorque
data.SrpTorque.sun_vector = [1;0;0];
Alr= 0.1*0.9*0.1-0.01*0.9*0.1;
data.SrpTorque.Alr = [0,0,0;
                     0,0,Alr;
                     0,Alr,0];
                 
%Configure MagRollControl
data.MagRollControl.A = [0.2,0.2,0.2];
data.MagRollControl.Y = 2012;
data.MagRollControl.k = 2e-4;
data.MagRollControl.wind_model=@WindHWM93;
data.MagRollControl.verb=1;

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

