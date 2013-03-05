% DDsatCo is a test script for OrbitProp
%
% In this case this example shows the orbit propagator in it's full
% potential using several models and simulating DDsat with the panels
% counter-rotating.
%
% In this case when the roll exceeds 1 deg the panels change orientation.
% This is done by stopping the integration (through an event) and changing
% the aerodynamic databases.
%
% The aerodynamic databases are in the AeroDB folder

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
data.sc_prop.I=[0.04,0,0;
                0,0.0177,0;
                0,0,0.0177]; %Inertia
            
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
data.models(end+1)={@MagDamping};
data.models(end+1)={@SrpTorque};
data.models(end+1)={@GravityGradient};
%data.models(end+1)={@MagRollControl};

%Configure AeroFDBCoWind
FP_rolldb=load(rolldb_FP);
FP_pitchdb=load(pitchdb_FP);
FP_yawdb=load(yawdb_FP);
data.AeroFDBCoWind.rolldb=FP_rolldb.F_roll;
data.AeroFDBCoWind.pitchdb=FP_pitchdb.F_pitch;
data.AeroFDBCoWind.yawdb=FP_yawdb.F_yaw;
data.AeroFDBCoWind.atmos_model=@nrlmsise00;
data.AeroFDBCoWind.wind_model=@WindHWM93;
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
data.AeroTDBCoWind.wind_model=@WindHWM93;
data.AeroTDBCoWind.CfA0l0=[TP_rolldb.T_roll(1,2);TP_rolldb.T_roll(1,3);TP_rolldb.T_roll(1,4)];

%Configure WindHWM93
data.WindHWM93.day0 = 10; %Initial day of the year
data.WindHWM93.solar_activity = 'Medium'; %Level of solar activity

%Configure MagDamping
data.MagDamping.c = 2e-4;
data.MagDamping.A = [0.2,0.2,0.2];
data.MagDamping.Y = 2012;
data.MagDamping.w = [0,deg2rad(-360/90/60),0];
data.MagDamping.verb = 0; 

%Configure SrpTorque
data.SrpTorque.sun_vector = [1;0;0];
Alr= 0.1*0.9*0.1-0.01*0.9*0.1;
data.SrpTorque.Alr = [0,0,0;
                     0,0,Alr;
                     0,Alr,0];
                 
%Configure MagRollControl
data.MagRollControl.A = [0.2,0.2,0.2];
data.MagRollControl.Y = 2012;
data.MagRollControl.k = 0;2e-4;
data.MagRollControl.wind_model=@WindHWM93;
data.MagRollControl.verb = 0;

%Events
data.event.function=@EvAtt;
data.EvAtt.wind_model=@WindHWM93;
data.EvAtt.r1=deg2rad(10);
data.EvAtt.r2=deg2rad(90);
data.EvAtt.r3=deg2rad(90);

%Configure CounterRotRollControl
data.CounterRotRollControl.wind_model=@WindHWM93;
data.CounterRotRollControl.FModel='AeroFDBCoWind';
data.CounterRotRollControl.TModel='AeroTDBCoWind';
data.CounterRotRollControl.Froll='rolldb';
data.CounterRotRollControl.Fpitch='pitchdb';
data.CounterRotRollControl.Fyaw='yawdb';
data.CounterRotRollControl.Troll='rolldb';
data.CounterRotRollControl.Tpitch='pitchdb';
data.CounterRotRollControl.Tyaw='yawdb';
FM_rolldb=load(rolldb_FM);
FM_pitchdb=load(pitchdb_FM);
FM_yawdb=load(yawdb_FM);
TM_rolldb=load(rolldb_TM);
TM_pitchdb=load(pitchdb_TM);
TM_yawdb=load(yawdb_TM);
data.CounterRotRollControl.FrollDbP=FP_rolldb.F_roll;
data.CounterRotRollControl.FpitchDbP=FP_pitchdb.F_pitch;
data.CounterRotRollControl.FyawDbP=FP_yawdb.F_yaw;
data.CounterRotRollControl.TrollDbP=TP_rolldb.T_roll;
data.CounterRotRollControl.TpitchDbP=TP_pitchdb.T_pitch;
data.CounterRotRollControl.TyawDbP=TP_yawdb.T_yaw;
data.CounterRotRollControl.FrollDbM=FM_rolldb.F_roll;
data.CounterRotRollControl.FpitchDbM=FM_pitchdb.F_pitch;
data.CounterRotRollControl.FyawDbM=FM_yawdb.F_yaw;
data.CounterRotRollControl.TrollDbM=TM_rolldb.T_roll;
data.CounterRotRollControl.TpitchDbM=TM_pitchdb.T_pitch;
data.CounterRotRollControl.TyawDbM=TM_yawdb.T_yaw;
data.CounterRotRollControl.CfA0='CfA0';
data.CounterRotRollControl.CfA0l0='CfA0l0';

%--- Integrate ---%
t_aux=tf(1);
t_f=[];
x_f=[];
tf_end=tf(end);
while t_aux(end)<tf(end)
    %Integrate
    [t_aux,x_aux] = OrbitProp(x0,tf(1),tf(2:end),data);
    
    %--- Attitude in flow reference frame ---%
    t=t_aux(end);
    t=t(:);
    y=x_aux(end,:);
    y=y(:);
    %Flow is velocity - atmosphere rotatio - wind
    wind = data.EvAtt.wind_model(t,y,data); %Wind
    we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
    corot = cross([0;0;we],y(1:3)); %Atmosphere rotation
    flow = y(4:6) - corot - wind;
    %Flow reference frame
    xa = flow/norm(flow); %Velocity (roll axis)
    zaux = -y(1:3)/norm(y(1:3)); %Nadir pointing
    ya = cross(zaux,xa); % Velocity anti-paralell (pitch axis)
    za = cross(xa,ya); % Near nadir (yaw axis)
    
    %Rotation matrix from ECEF to flow ref
    Rfe=[xa';ya';za'];
    
    %Rot from ECEF to body
    Rbe=quat2dcm(y(7:10)');
    
    %Rot from flow to body
    Rbf = Rbe*inv(Rfe);
    %Get angles
    [r1, r2, r3] = dcm2angle(Rbf,'XYZ');
    
    %Change the aerodatabases.
    if r1>0
        %Roll higher than the treshold
        %Change databeses
        %Force
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Froll,data.CounterRotRollControl.FrollDbM);
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Fpitch,data.CounterRotRollControl.FpitchDbM);
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Fyaw,data.CounterRotRollControl.FyawDbM);
        CfA0=[data.CounterRotRollControl.FrollDbM(1,2);data.CounterRotRollControl.FrollDbM(1,3);data.CounterRotRollControl.FrollDbM(1,4)];
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.CfA0,CfA0);
        %Torque
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Troll,data.CounterRotRollControl.TrollDbM);
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Tpitch,data.CounterRotRollControl.TpitchDbM);
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Tyaw,data.CounterRotRollControl.TyawDbM);
        CfA0l0=[data.CounterRotRollControl.TrollDbM(1,2);data.CounterRotRollControl.TrollDbM(1,3);data.CounterRotRollControl.TrollDbM(1,4)];
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.CfA0l0,CfA0l0);
        
    elseif r1<0
        %Roll lower than the treshold
        %Change databeses
        %Force
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Froll,data.CounterRotRollControl.FrollDbP);
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Fpitch,data.CounterRotRollControl.FpitchDbP);
        data=setfield(data,data.CounterRotRollControl.FModel,data.CounterRotRollControl.Fyaw,data.CounterRotRollControl.FyawDbP);
        CfA0=[data.CounterRotRollControl.FrollDbP(1,2);data.CounterRotRollControl.FrollDbP(1,3);data.CounterRotRollControl.FrollDbP(1,4)];
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.CfA0,CfA0);
        %Torque
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Troll,data.CounterRotRollControl.TrollDbP);
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Tpitch,data.CounterRotRollControl.TpitchDbP);
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.Tyaw,data.CounterRotRollControl.TyawDbP);
        CfA0l0=[data.CounterRotRollControl.TrollDbP(1,2);data.CounterRotRollControl.TrollDbP(1,3);data.CounterRotRollControl.TrollDbP(1,4)];
        data=setfield(data,data.CounterRotRollControl.TModel,data.CounterRotRollControl.CfA0l0,CfA0l0);
    end
        
    %Save time stamps and state vector
    t_f=[t_f;t_aux];
    x_f=[x_f;x_aux];
    
    %Initial conditions
    x0=x_f(end,:);
    tf=t_f(end):1:tf_end;
    
    fprintf('Changing panel deflections!\n');
    fprintf('Time between changes: %g min\n',(t_aux(end)-t_aux(1))/60);
    fprintf('Current time: %g min\n',t_aux(end)/60);
    
    %Change event limit
    data.EvAtt.r1=-data.EvAtt.r1;
end
t=t_f;
x=x_f;

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
 [A,Torque(end+1,1:3),E] = MagDamping(t(i),x(i,:)',data);
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
    mfield_ECEF = igrf11magm(h, lat, lon, data.MagDamping.Y)/1e9;
    %Change magnetic field to body axes.
    DCM = quat2dcm([x(i,10),x(i,7:9)]);
    mfield = DCM*mfield_ECEF';
    
    %Preliminary actuation levels m_p
    m_p=(cross(mfield,Torque(i,:))/norm(mfield)^2)';
    %Check different scalar values
    k_max=(data.MagDamping.A'-m_p)./mfield; %Scalar value for maximum actuation level on one axis
    k_min=(-data.MagDamping.A'-m_p)./mfield; %Scalar value for minus maximum actuation level on one axis
    k_0=-m_p./mfield; %Scalar value for 0 actuation level on one axis
    
    %Scan through the different scalar values
    err=1e-3/100;
    m_v=[];
    p_v=[];
    for k=[k_max',k_min',k_0']
        m=m_p+k*mfield;
        
        if abs(m(1))>data.MagDamping.A(1)*(1+err)
            %Solution is not real
            continue
        end
        if abs(m(2))>data.MagDamping.A(2)*(1+err)
            %Solution is not real
            continue
        end
        if abs(m(3))>data.MagDamping.A(3)*(1+err)
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

