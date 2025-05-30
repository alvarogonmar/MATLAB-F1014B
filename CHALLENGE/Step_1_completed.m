%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session 1 | Step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rodrigo Gamboa & Francisco Montes | May 2024

% NO USAR ! ! !

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
                           %this for: almacenar los componentes del campo magnético generado por cada dl en cada punto del espacio 3D

for I = 1:Lx               %Open a for loop using "I" as index, going from 1 to Lx.  
    for J = 1:Ly           %Same but for Ly, use a "J" index.
        for K = 1:Lz       %Same but for Lz, use a "K" index.   

            for L = 1:(nl*N)   %Open a 4th for loop, usin an "L" index, running from 1 to nl*N (Why? porque es el número total de elementos de corriente)
            
                % r = R - r'
                rx = x(I) - Px(L);     %Write here the rx component
                ry = y(J) - Py(L);     %Same but for the ry component
                rz = z(K) - Pz(L);     %Same but for the rz component

                r = sqrt(rx^2 + ry^2 + rz^2);   %Get the magnitude

                r3 = r^3;              %Declate an r3 variable
                                       %equal to the third power of r
                
                % Cross product dl x r = (dy*rz - dz*ry)i + (dz*rx - dx*rz)j + (dx*ry - dy*rx)k
                dBx(I,J,K) = dBx(I,J,K) + km * (dy(L)*rz - dz(L)*ry)/r3;    %Write here the dBx(I,J,K) component
                dBy(I,J,K) = dBy(I,J,K) + km * (dz(L)*rx - dx(L)*rz)/r3;    %Same but for the dBy(I,J,K) component  
                dBz(I,J,K) = dBz(I,J,K) + km * (dx(L)*ry - dy(L)*rx)/r3;    %Same but for the dBz(I,J,K) component

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

