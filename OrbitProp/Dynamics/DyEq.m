function [ydot] = DyEq(t,y,data)
% DyEq computes the derivates of the state vector.
%
% For derivates of position and attitude quaternion uses linear and angular
% velocity.
%
% For derivates of linear and angular velocity uses the specified models.
% Ecach model must output two 3x1 vectors the first one containing the
% linear acceleration in m/s2 and the second the torque in Nm.
% The models will be passed the following parameters (time, state_vector,
% data).
%
% It is considered good practice to store the model dependant parameters
% within a structure in data named after the model name (i.e
% data.model_name.parameter).

%--- Copyright notice ---%
% Copyright 2012 Cranfield University
% Written by Josep Virgili
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

%--- Velocity ---%
%Linear velocity is the derivate of position
ydot(1:3) = y(4:6);

%Compute angular acceleration matrix
w = y(11:13);
Om=[0 w(3) -w(2) w(1);
    -w(3) 0 w(1) w(2);
    w(2) -w(1) 0 w(3);
    -w(1) -w(2) -w(3) 0];

%Derivate of the quaternion
qdot = 1/2*Om*[y(8:10);y(7)];
%The formula above is for scalar part last. Sort it.
ydot(7)=qdot(4);
ydot(8:10) = qdot(1:3);


%--- Acceleration ---%
%Initialize variables
A=0;
T=0;
E=0;

%Call models
for i=1:size(data.models,1)
    %Make the call
    [Aaux,Taux,Eaux] = data.models{i}(t,y,data);
     
    %Add the model results to the total 
    A=A+Aaux;
    T=T+Taux;
    E=E+Eaux;
end

%Linear acceleration
ydot(4:6) = A;

%Angular accelaration
ydot(11:13)=data.sc_prop.I\(T-cross(y(11:13),data.sc_prop.I*y(11:13)));

%Other state variables and format output
ydot=[ydot(:);E(:)];

end