% This mfile is to calculate momentum initial matrix for DDsat
%
% %--- Copyright notice ---%
% Copyright 2012-2013 Cranfield University
% Written by Daniel Zhou Hao
%
% This file is part of the massive attitude control simulation for DDsat
%
% Mess budget: Unit in kg
%
m_total = 2.380;
m_all_aerofoil = 0.3356;
%
m_each_aerofoil = m_all_aerofoil /4;
m_cube_body = m_total - m_all_aerofoil;
% Structure Size: Unit in m
%
% 2U Body
 d_body = 0.1; % Depth
 w_body = 0.1; % Width
 h_body = 0.225; % Height
% Aerofoil 
h_aerofoil = 0.180*3 ; % Height
w_aerofoil = 0.065 ; % Width
s_aerofoil = h_aerofoil * w_aerofoil;
r_xx= h_aerofoil / 2 + w_body / 2; % Inertia Shift from aerofoil centre to body contre about xx axis
r_yy =  h_body/2 - (w_aerofoil / 2);
r_zz = r_yy;
 % Apply Momentum Intertia Paraller Aixs Therom
Ixx = (1/12) * m_cube_body* ((d_body)^2 + (w_body)^2) + 4 * m_each_aerofoil * (w_aerofoil *( h_aerofoil ^ 3 ) /12 + s_aerofoil * (r_xx^2));

% Apply Momentum Perpendicular Axis Theorem;
% Iz = Ix + Iy, where Ix and Iy is in on the same plane as the rotational
% plane: Therefore, for the aerofoil rotates along with a prependicular
% axis, Izz = (h_aerofoil *( w_aerofoil ^ 3 ) /12 + s_aerofoil * (r_yy^2))
% +  (w_aerofoil *( h_aerofoil ^ 3 ) /12 + s_aerofoil * (r_xx^2))
Iyy = (1/12)  * m_cube_body* ((h_body)^2 + (d_body)^2)  + m_each_aerofoil *(2*  ((h_aerofoil *( w_aerofoil ^ 3 ) /12 + s_aerofoil * (r_yy^2))  +  (w_aerofoil *( h_aerofoil ^ 3 ) /12 + s_aerofoil * (r_xx^2))));
Izz = Iyy;
I = [Ixx,Iyy,Izz]

% For PDR: 
% I =
%     0.0049    0.0111    0.0111
