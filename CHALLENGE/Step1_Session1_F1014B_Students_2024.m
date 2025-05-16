%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session 1 | Step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rodrigo Gamboa & Francisco Montes | May 2024


%-OBJECTIVE-%

%Create and visualize a discretized loop made of N arrows, 
%simbolizing electric current differentials, so that the whole 
%configuration will approximate a circular wire (coil), with an 
%associated constant electric current I. The more currents the better.

%Recall commenting and completing all necessary code from now on,
%as this will impact the quality and grade of the associate Deliverable.


%-----------------------------------------------------------------------
%Code starts here
%Complete the code in the left, using the comments in the right.
%-----------------------------------------------------------------------

clc;          %Clear all variables from the workspace
clear;        %Clear command window
close all;    %Clear all elements on figures

nl = 1;              %Declare a variable "nl" that will
                     %represent the number of wire loops in your
                     %simulation.One for a configuration called 
                     % a ______________? %Two or more for a configuration 
                     % referred to as a __________________? 

ds = 0.1;             % Declare a variable "ds", for lenght differential.
                      % Set it equal to 0.1 as a starting point.

x = -5:ds:5; y = x;  z = x;                      %Define x, y and z vectors, all equal, 
                                                 %ranging from -5 to 5, in steps of ds.
                                                 %This will define our 3D working space.
                                    
                                                 %When correctly done, each vector will
                                                 %appear in the Workspace as a 1x101
                                                 %double array.

Lx = length(x); Ly = length(y); Lz = length(z);                                           %Using the length(A) function, 
                                                                                          %where A can be a vector or array, 
                                                                                          %declare Lx, Ly and Lz
                                                                                          %variables, setting them
                                                                                          %equal to the number of elements 
                                                                                          %in the x, y and z vectors (i.e. the 
                                                                                          %vector lenghts).  

rw = 0.2;                                         %Using a variable labelled rw, declare the 
                                                 %wire thickness, use a value of 0.2.

I = 300;                                          %Using a variable I, declare the 
                                                 %associate electric current
                                                 %(use value of 300 Amperes).

%Uncomment the following code.

mo=4*pi*1e-7;                                   %What is this? % Permeabilidad del vacío [T·m/A] constante
km=mo*I/(4*pi);                                 %And this?     % Constante magnética simplificada

N = 20;                       % (flechas) Declare a variable N, that will define the # of points per loop.
R = 1.5;                       %Declare a variable R, as the radious if the wire, use 1.5 as default.                   
sz = 1;                        %Declare a variable sz, setting it equal to 1. This will represent the 
                              %loop step size in the z axis direction.

s = 1;                         %Declare an s variable (equal to 1), representing the loop "number".

dtheta = 2*pi/N;              %Define the "differential angle step" (dtheta) with respect to the center 
                              %of the circular wire (coil), consider that it has to cover an
                              %angular displacement of 2*pi. Think about it!

dl = R * dtheta;              %Define the "dl" length differential, coupled with the associated dtheta.   
                              %Consider that dtheta is in radians.

ang = 0:dtheta:(2*pi - dtheta);    %Define a vector called "ang" (for angle), 
                                   %ranging from 0 to 2*pi-dtheta, in steps of dtheta.
                                   %This will store angle values to complete a cycle 
                                   % with steps dtheta, around the coil. Why from 
                                   % 0 to 2*pi-dtheta:_______________________________? 

%Start here a for loop, looping a dummie variable
%going from 1 to nl (the number of loops in the system).

% Posiciones x
Px = zeros(1, N*nl);
Py = zeros(1, N*nl);
Pz = zeros(1, N*nl);
% Componentes x,y,z de la flecha
dx = zeros(1, N*nl);
dy = zeros(1, N*nl);
dz = zeros(1, N*nl);

for i = 1:nl

    %Px(s:s+N-1)=R*cos(ang);           %Comment
    %Py(s:s+N-1)=R*sin(ang);           %Comment
    %Pz(s:s+N-1)=-nl/2*sz+(i-1)*sz;    %Comment

    %dx(s:s+N-1)=-Py(s:s+N-1)*dtheta;  %Comment
    %dy(s:s+N-1)=Px(s:s+N-1)*dtheta;   %Comment

    %s=s+N;                            %Comment

    Px(s:s+N-1)=R*cos(ang);           
    Py(s:s+N-1)=R*sin(ang);           
    Pz(s:s+N-1)=-nl/2*sz+(i-1)*sz;    

    dx(s:s+N-1)=-Py(s:s+N-1)*dtheta;  % dx ~ -sin → dirección tangente al círculo
    dy(s:s+N-1)=Px(s:s+N-1)*dtheta;   % dy ~ cos

    s = s + N; % Se avanza al siguiente bloque para el siguiente loop 

end    

dz(1:N*nl)=0;                          %Comment        

%Call figure 1
figure(1)

%Use quiver3 function, to visualize the wire of currents.
%Recall that quiver3(X,Y,Z,U,V,W) plots arrows with directional components 
%U, V, and W at the Cartesian coordinates specified by X, Y, and Z.

quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2)  %What are the other options for?

%Uncomment for setting a proper viewing perspective.
view(-34,33)
