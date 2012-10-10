% Drag_J4 is a test script for OrbitProp
%
% This is script uses the a gravity model with up to the J4 harmonics. Also
% it includes drag using a simple ballistic coefficient
%
% The integration stops when the spacecarft re-enters (assumed to be at 100
% km).

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
h = 250*1e3; %Initial altitude in m
i = 70; %Inclination
v0 = sqrt(mu/(Re+h)); %Initial velocity
x0 = [Re+h,0,0,0,v0*cosd(i),v0*sind(i)]; %Initial vector state
tf = inf; %See event for integration stopping condition

%Models
data.models={@GravJ4,@ConstantBCCoWind};

%Configure ConstantBC model
data.ConstantBCCoWind.BC=2/(0.1*0.1*2.5); %Ballisic coefficient
data.ConstantBCCoWind.atmos_model=@nrlmsise00;
data.ConstantBCCoWind.wind_model=@NoWind;

%Configure nrlmsise00
data.nrlmsise00.solar_activity = 'Medium'; %Level of solar activity
data.nrlmsise00.day0 = 10; %Initial day of the year

%Event
data.event.function = @EvAlt;
data.event.h = 100e3;

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
plot(t,lla(:,3)/1e3)
xlabel('Time [min]')
ylabel('Altitude [km]')