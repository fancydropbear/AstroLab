% SphericalGrav is a test script for OrbitProp
%
% This is script shows the differences between the different gravity
% models.

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

%Initial conditions
h = 200*1e3; %Initial altitude in m
i = 70; %Inclination
v0 = sqrt(mu/(Re+h)); %Initial velocity
x0 = [Re+h,0,0,0,v0*cosd(i),v0*sind(i)]; %Initial vector state
tf = 1 * 90 * 60; %Integration time in s

%--- Iterate for different Gravity models---%
%Spherical gravity
data.models={@GravSpherical};
%Integrate
[t,x] = OrbitProp(x0,0,tf,data);
%Compute altitude
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
%Save relevant variables
hsphe=lla(:,3)/1e3;
tsphe=t;

%J4 gravity
data.models={@GravJ4};
%Integrate
[t,x] = OrbitProp(x0,0,tf,data);
%Compute altitude
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
%Save relevant variables
hJ4=lla(:,3)/1e3;
tJ4=t;

%EGM96 gravity
data.models={@GravEGM96};
%Integrate
[t,x] = OrbitProp(x0,0,tf,data);
%Compute altitude
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
%Save relevant variables
hEGM96=lla(:,3)/1e3;
tEGM96=t;

%EGM2008 gravity
data.models={@GravEGM2008};
%Integrate
[t,x] = OrbitProp(x0,0,tf,data);
%Compute altitude
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
%Save relevant variables
hEGM2008=lla(:,3)/1e3;
tEGM2008=t;

%--- Plot Altitudes ---%
figure(1)
plot(tsphe/60,hsphe,tJ4/60,hJ4,tEGM96/60,hEGM96,tEGM2008/60,hEGM2008)
legend('Spherical','J4','EGM96','EGM2008')
xlabel('Time [min]')
ylabel('Altitude [km]')
