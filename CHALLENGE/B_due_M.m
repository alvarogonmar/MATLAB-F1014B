%-------------------------------%
%   Challenge Session No. 4
%-------------------------------%
%Rodrigo Gamboa & Francisco Montes | May 2024


%------------Part 3---------------% 
%FUNCTION DEFINITION | DECLARATION
%---------------------------------%

%Change the name of this tab to the name of the referred or associated
%function (the one that got you here in the first place). 
%e.g. If the function called in the main code is f_meg(x,y), this
%tab should be called then f_meg.m, for Matlab to locate it. Be sure to 
%save it in the same directory as your main code. 

%Definition or declaration of a function called
%[x,y,phiB2,Bz]=B_due_M(zm(cc+1),mag,Rring);
%that takes (input) _____, ____ and _____; and returns (output)
%_____, ______, ______ and _____.

%Declare the function here, start by checking reserved word/code for
%MatLab, to declare a function.
function [x, y, phiB, Bz] = B_due_M(z, mag, Rring)

%[Reserved word for definig a function] <function call>


mo=4*pi*1e-7;           % Constante de permeabilidad magnética del vacío (μ₀), en H/m.
ds=0.005;               %cada cuadrito de la cuadrícula que cubre el área de la espira mide 0.005 metros de lado. 
                        % Sirve para dividir el área total en pedacitos pequeños y sumar el campo magnético en cada uno.

x = -Rring:ds:Rring;    % Vector de coordenadas x que cubre el área de la espira                 
y = -Rring:ds:Rring;    % Vector de coordenadas y que cubre el área de la espira
               
Lx = length(x);         % Número de puntos en el eje x
Ly = length(y);         % Número de puntos en el eje y

Bz = zeros(Lx, Ly);     % Inicializa una matriz Bz de ceros para guardar el campo en cada punto (x,y)
phiB = 0;               % Inicializa el flujo magnético total en cero


% Doble ciclo para recorrer todos los puntos de la cuadrícula en el plano xy
for i=1:Lx          
    for j=1:Ly

        r=sqrt(x(i)^2+y(j)^2);       % Calcula la distancia desde el centro (radio) hasta el punto (x(i), y(j))
                                     % Esto es necesario para saber si el punto está dentro del área circular de la espira.

        if r<Rring             % Si el punto está dentro del radio de la espira, se calcula el campo Bz ahí


            % Fórmula del campo magnético vertical Bz de un dipolo en coordenadas
            Bz(i,j) = mo/(4*pi) * mag * (3*z^2 - (x(i)^2 + y(j)^2 + z^2)) / ((x(i)^2 + y(j)^2 + z^2)^(5/2));

            % "STAGE 3: Calculation of induced fem by a magnet (dipole) free falling through a coil"

            phiB= phiB + Bz(i, j) * ds^2;  % phiB acumula el flujo magnético total (ΦB) que atraviesa la espira
                                           % usando la suma de Bz * área en cada punto

        end
    end
end