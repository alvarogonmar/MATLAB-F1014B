%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session 1 | Step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rodrigo Gamboa & Francisco Montes | May 2024

%-OBJECTIVE-%
% Create and visualize a discretized loop made of N arrows, 
% simbolizing electric current differentials, so that the whole 
% configuration will approximate a circular wire (coil), with an 
% associated constant electric current I. The more currents the better.

%-----------------------------------------------------------------------
% Code starts here
%-----------------------------------------------------------------------

clc;                                  % Clear all variables from the workspace
clear;                                % Clear command window
close all;                            % Clear all elements on figures

nl = 1;                               % Number of wire loops. One = loop (espira),
                                      % Two or more = solenoid (bobina).

ds = 0.1;                             % Length differential

x = -5:ds:5; y = x; z = x;            % Define x, y, z vectors from -5 to 5 in steps of ds

Lx = length(x);                      % Number of elements in x
Ly = length(y);                      % Number of elements in y
Lz = length(z);                      % Number of elements in z

rw = 0.2;                             % Wire thickness
I = 300;                              % Electric current in Amperes

mo = 4*pi*1e-7;                       % Permeability of free space (H/m)
km = mo*I/(4*pi);                    % Constant for Biot–Savart Law (units: T·m/A)

N = 100;                              % Number of points per loop
R = 1.5;                              % Radius of the wire
sz = 1;                               % Step size in z-direction

s = 1;                                % Starting index for points in vector

dtheta = 2*pi/N;                     % Angle step for a full circle (radians)
dl = R * dtheta;                     % Arc length differential

ang = 0:dtheta:(2*pi - dtheta);      % Angle values to complete a cycle
                                     % 2*pi - dtheta avoids repeating the start point

% Initialize vectors for positions and differentials
Px = zeros(1, N*nl);
Py = zeros(1, N*nl);
Pz = zeros(1, N*nl);
dx = zeros(1, N*nl);
dy = zeros(1, N*nl);
dz = zeros(1, N*nl);

for i = 1:nl                         % Loop over the number of wire loops
    Px(s:s+N-1) = R * cos(ang);                     % x-coordinates of the loop
    Py(s:s+N-1) = R * sin(ang);                     % y-coordinates of the loop
    Pz(s:s+N-1) = -nl/2*sz + (i-1)*sz;              % z-plane of the loop
    
    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;            % x-component of dl (tangent vector)
    dy(s:s+N-1) = Px(s:s+N-1) * dtheta;             % y-component of dl (tangent vector)

    s = s + N;                          % Update index for next loop
end

dz(1:N*nl) = 0;                        % No z-component of current (in-plane)

figure(1)                             % Open figure window

% Plot the current vectors using quiver3
quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2)  
% 0.5 scales arrows, '-r' makes them red, LineWidth=2 makes them thicker

view(-34, 33)                         % Set viewing perspective
axis equal                            % Keep proportions in 3D
xlabel('X'); ylabel('Y'); zlabel('Z'); % Axis labels
title('Circular Loop of Electric Current')
grid on
