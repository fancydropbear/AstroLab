% SphericalGrav is a test script for OrbitProp
%
% This is script propagates an orbit using only an spherical gravity model.

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
tf = 10 * 90 * 60; %Integration time in s

%Models
data.models={@GravSpherical};

%--- Call function ---%
[t,x] = OrbitProp(x0,0,tf,data);

%--- Post process ---%
%Compute altitude lat and long
we = 7.292115e-5; %Earth Angular velocity in rad/s [source: SMAD 3rd edition]
Req=6378136.49;   %Equatorial Earth radius m [source: SMAD 3rd edition]
f=1/298.256;      %Flattening factor [source: SMAD 3rd edition]

%Adjust Longitude
lla=ecef2lla(x(:,1:3),f,Req); %Compute taking assuming Earth as an ellipsoid
lla(:,2) = lla(:,2)- t*we*180/pi;

%Plot ground track
figure(1)
geoshow('landareas.shp', 'FaceColor', [0.5 1.0 0.5]);
axis tight
geoshow(lla(:,1),lla(:,2));

%Plot Altitude
figure(2)
plot(t/60,lla(:,3)/1e3)
xlabel('Time [min]')
ylabel('Altitude [km]')