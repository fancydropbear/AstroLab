% Statics for Bdot Detumbling

% Statics for the Detumbling
% This function is dedicated to get the (detumbling time statics), (power
% statics)
% Before running this static file, remove the 'clear function' from th
% Detumbling.m

clear all;
clc;


for counter = 1:30
    
    % Find where the time for angular rate getting to small enough
    DDsatDetumblingBDotError;
    
    for counter2 = 1: length(t)
        if norm(rad2deg(x(counter2,11:13))) < 0.5
            detumb_t(counter) = t(counter2);
            avg_power(counter) = mean(PT(1:counter2));
            break;
        end
    end
end

mean_t = mean(detumb_t(:))/60
std_t = std(detumb_t(:))/60
mean_avg_power = mean(avg_power(:))
std_avg_power = std(avg_power(:))
