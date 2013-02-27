% DDsatDetumbling is a test script for OrbitProp
%
% In this case this example shows a detumbling example for DDsat using
% magnetic torquers

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

tf = 0:1:10 * 60; %Integration time
rx=deg2rad(-90);
ry=deg2rad(-90);
rz=deg2rad(-atand((v0*sind(i))/(v0*cosd(i)-Vco)));
q0=angle2quat(rx,ry,rz,'XYZ'); %Initial attitude quaterion

%Initial angular rate (make it random!)
w1=20*(rand()-0.5);
w2=20*(rand()-0.5);
w3=20*(rand()-0.5);
w0=deg2rad([w1,w2,w3]); %Format initial angular rate

%Format the initial 
x0=[x0,q0,w0];

%Spacecarft propetries
data.sc_prop.I=[0.0042,0,0;
                0,0.0104,0;
                0,0,0.0104]; %Inertia
            
%Models
data.models={@GravJ4};
data.models(end+1)={@MagDamping};

%Configure MagDamping
data.MagDamping.c = 2e-3;
data.MagDamping.A = [0.2,0.2,0.2];
data.MagDamping.Y = 2012;
data.MagDamping.w = [0,deg2rad(-360/90/60),0];
data.MagDamping.verb=1;

%--- Integrate ---%
[t,x] = OrbitProp(x0,tf(1),tf(2:end),data);

%--- POST-PROCESS ---%

%--- Rate plot ---%
figure
plot(t/60,rad2deg(x(:,11)),t/60,rad2deg(x(:,12)),t/60,rad2deg(x(:,13)));
xlabel('Time [min]')
ylabel('Angular rates [deg/s]')
legend('Roll','Pitch','Yaw')
title('Angular rates')

