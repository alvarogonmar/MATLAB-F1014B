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

for I = 1:nl

    %Px(s:s+N-1)=R*cos(ang);           %Comment
    %Py(s:s+N-1)=R*sin(ang);           %Comment
    %Pz(s:s+N-1)=-nl/2*sz+(i-1)*sz;    %Comment

    %dx(s:s+N-1)=-Py(s:s+N-1)*dtheta;  %Comment
    %dy(s:s+N-1)=Px(s:s+N-1)*dtheta;   %Comment

    %s=s+N;                            %Comment

    Px(s:s+N-1)=R*cos(ang);           
    Py(s:s+N-1)=R*sin(ang);           
    Pz(s:s+N-1)=-nl/2*sz+(I-1)*sz;    

    dx(s:s+N-1)=-Py(s:s+N-1)*dtheta;  % dx ~ -sin → dirección tangente al círculo
    dy(s:s+N-1)=Px(s:s+N-1)*dtheta;   % dy ~ cos

    s = s + N; % Se avanza al siguiente bloque para el siguiente loop 

end    

dz = zeros(1, N*nl);                          %Comment        

%Call figure 1
figure(1)

%Use quiver3 function, to visualize the wire of currents.
%Recall that quiver3(X,Y,Z,U,V,W) plots arrows with directional components 
%U, V, and W at the Cartesian coordinates specified by X, Y, and Z.

quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2)  %What are the other options for?

%Uncomment for setting a proper viewing perspective.
view(-34,33)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic-field Calculation & Visualization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic-field Calculation & Visualization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Using the results obtained from the Biot-Savart law, get the B-field 
%at each point in space, around the current loop.

dBx = zeros(Lx, Ly, Lz);   %Create three, 3D equal arrays, dBx,
dBy = zeros(Lx, Ly, Lz);   %dBy & dBz, each spanning from 1
dBz = zeros(Lx, Ly, Lz);   %to Lx, 1 to Ly and 1 to Lz, and 
                           %initialize all values to zero. You should 
                           %obtain in the workspace, three 101x101x101 
                           % objects (dBx, dBy, dBz). 
                           
                           %What do you think we are to use
                           %this for: almascenar los componentes del campo magnético generado por cada dl en cada punto del espacio 3D

for I = 1:Lx               %Open a for loop using "I" as index, going from 1 to Lx.  
    for J = 1:Ly           %Same but for Ly, use a "J" index.
        for K = 1:Lz       %Same but for Lz, use a "K" index.   

            for L = 1:(nl*N)   %Open a 4th for loop, usin an "L" index, running from 1 to nl*N (Why? porque es el número total de elementos de corriente)
            
                % r = R - r'
                rx = x(I) - Px(L);     %Write here the rx component
                ry = y(J) - Py(L);     %Same but for the ry component
                rz = z(K) - Pz(L);     %Same but for the rz component

                r = sqrt(rx^2 + ry^2 + rz^2+rw^2);   %Get the magnitude

                r3 = r^3;              %Declate an r3 variable
                                       %equal to the third power of r
                
                % Cross product dl x r = (dy*rz - dz*ry)i + (dz*rx - dx*rz)j + (dx*ry - dy*rx)k
                % dBx(I,J,K) = dBx(I,J,K) + km * (dy(L)*rz - dz(L)*ry)/r3;    %Write here the dBx(I,J,K) component
                % dBy(I,J,K) = dBy(I,J,K) + km * (dz(L)*rx - dx(L)*rz)/r3;    %Same but for the dBy(I,J,K) component  
                % dBz(I,J,K) = dBz(I,J,K) + km * (dx(L)*ry - dy(L)*rx)/r3;    %Same but for the dBz(I,J,K) component
                dBx (I,J,K) = dBx(I, J, K) + km * dy(L) * rz / r3;
                dBy (I,J,K) = dBy(I, J, K) + km * dx(L) * rz / r3;
                dBz (I,J,K) = dBz(I, J, K) + km * (dx(L) *ry - dy(L) * rx) / r3;

            end   %Close the 4th loop

        end    %Close the 3rd loop   
    end   %Close the 2nd loop
end    %Close the 1st loop

%------------------------------------------------------------
%Uncomment all the following, and explain (comment) each line
%------------------------------------------------------------

Bmag = sqrt(dBx.^2 + dBy.^2 + dBz.^2);          % Magnitud del campo B en cada punto 3D

centery = round(Ly/2);                          % Selecciona la mitad del plano Y para hacer el corte XZ
Bx_xz = squeeze(dBx(:,centery,:));              % Extrae los valores de Bx en ese plano (matriz 2D XZ)
Bz_xz = squeeze(dBz(:,centery,:));              % Extrae los valores de Bz en ese plano (matriz 2D XZ)
Bxz = squeeze(Bmag(:,centery,:));               % Extrae la magnitud de B en ese plano

figure(2)                                       % Crea una nueva figura para visualizar el campo

hold on                                         % Mantiene la gráfica activa para múltiples elementos

pcolor(x, z, (Bxz').^(1/3)); shading interp; colormap jet; colorbar   
% Crea un mapa de color con la magnitud de B^(1/3) para mejor visualización;
% shading interp suaviza el gráfico; colormap jet aplica una escala de color;
% colorbar muestra la escala de colores

h1 = streamslice(x, z, Bx_xz', Bz_xz', 3);      % Dibuja líneas de corriente del campo magnético (streamlines)
set(h1, 'Color', [0.8 1 0.9]);                  % Cambia el color de las líneas a verde claro

xlabel 'x'                                      % Etiqueta eje x
ylabel 'z'                                      % Etiqueta eje z
title 'Magnetic field of a circular current'    % Título del gráfico
zz(length(ang))=0;
