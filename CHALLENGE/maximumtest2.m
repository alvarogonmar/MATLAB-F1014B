%---------------------------------------------------------%
%     Simulación de imán cayendo por una espira con FEM   %
%---------------------------------------------------------%
clc; clear; close all;  % Limpia la consola, las variables y las figuras

%------------------------ ESPIRA ------------------------%
nl = 1;        % Número de niveles (solo 1 espira)
N = 20;        % Número de segmentos en la espira (puntos en el círculo)
totalPoints = N * nl;  % Total de puntos a graficar

% Rango del espacio para visualización (malla)
ds = 0.1;              % Resolución del espacio
x = -5:ds:5;           % Coordenadas x
y = -5:ds:5;           % Coordenadas y
z = x;                 % Coordenadas z iguales a x (solo para generar el mismo tamaño)
Lx = length(x);        % Tamaño del eje x
Ly = length(y);        % Tamaño del eje y
Lz = length(z);        % Tamaño del eje z

% Parámetros físicos de la espira
rw = 0.2;              % Radio del alambre (no se usa directamente aquí)
I = 300;               % Corriente en la espira
mo = 4*pi*1e-7;        % Permeabilidad magnética del vacío
km = mo*I/(4*pi);      % Constante para campo magnético
R = 1.5;               % Radio de la espira
dtheta = 2*pi/N;       % Incremento angular para crear el círculo
ang = 0:dtheta:(2*pi - dtheta);  % Ángulos para cada punto en la espira

% Inicialización de vectores para posiciones y tangentes de la espira
Px = zeros(1, totalPoints);  % Coordenadas x de los puntos de la espira
Py = zeros(1, totalPoints);  % Coordenadas y
Pz = zeros(1, totalPoints);  % Coordenadas z (todas en plano z = 0)
dx = zeros(1, totalPoints);  % Componentes x de los vectores tangentes
dy = zeros(1, totalPoints);  % Componentes y
dz = zeros(1, totalPoints);  % Componentes z (cero porque está en el plano xy)

% Generación de los puntos de la espira
s = 1;  % Índice inicial
for I = 1:nl
    Px(s:s+N-1) = R * cos(ang);                 % Coordenadas x en círculo
    Py(s:s+N-1) = R * sin(ang);                 % Coordenadas y
    Pz(s:s+N-1) = 0;                            % z = 0 (plano xy)
    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;        % Derivada de x respecto a theta
    dy(s:s+N-1) = Px(s:s+N-1) * dtheta;         % Derivada de y respecto a theta
    s = s + N;                                  % Avanza al siguiente grupo
end

dz = zeros(1, totalPoints);  % Componentes z de los vectores tangentes = 0

%-------------------- PARÁMETROS DEL IMÁN --------------------%
mag = 500;       % Momento magnético del imán
Rring = 0.5;     % Radio del imán (suponiendo que también es circular)
zo = 0.1;        % Posición inicial del imán sobre la espira (eje z)
dt = 0.01;       % Intervalo de tiempo para la simulación
zring = 0;       % Altura de la espira (siempre en z = 0)

% Inicialización de vectores para simulación
t = zeros(1, 2000);      % Tiempo
zm = zeros(1, 2000);     % Posición del imán (z)
zm(1) = zo;              % Posición inicial del imán
vz = zeros(1, 2000);     % Velocidad del imán
fem = zeros(1, 2000);    % Fuerza electromotriz inducida
cc = 1;                  % Contador del tiempo

%------------------ SIMULACIÓN PRINCIPAL ------------------%
% Mientras el imán no haya pasado por completo la espira
while zm(cc) > 0.0172

    pause(0.001);  % Pausa breve para visualizar la animación

    %----------- ANIMACIÓN EN FIGURA 1 -----------%
    figure(1); clf                          % Limpia la figura
    quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2) % Dibuja la espira como vectores
    hold on
    scatter3(0, 0, zm(cc), 100, 'filled')  % Dibuja el imán como un punto grande en z
    axis([-2 2 -2 2 -0.2 0.2])              % Límites de la gráfica
    xlabel('x'); ylabel('y'); zlabel('z');
    title(['Posición del imán (z = ' num2str(zm(cc), '%.3f') ' m)'])  % Muestra posición actual
    view(3)                                % Vista en 3D
    drawnow                                % Actualiza la figura

    %----------- CÁLCULO DE LA FEM INDUCIDA -----------%
    [xg, yg, phiB1, Bz] = B_due_M(zm(cc), mag, Rring);     % Flujo magnético en t
    zm(cc+1) = zm(cc) + vz(cc)*dt - 0.5*9.81*dt^2;         % Nueva posición (fórmula de física)
    vz(cc+1) = vz(cc) - 9.81*dt;                           % Nueva velocidad (aceleración gravitacional)
    [~, ~, phiB2, ~] = B_due_M(zm(cc+1), mag, Rring);      % Flujo magnético en t+1
    fem(cc) = - (phiB2 - phiB1) / dt;                      % Ley de Faraday: FEM inducida

    %----------- GRAFICACIÓN COMENTADA (opcional) -----------%
    % Las siguientes líneas están comentadas, pero servirían para mostrar
    % gráficas adicionales como la FEM, altura, o el campo magnético en tiempo real

    % subplot(2,2,1)  → FEM en el tiempo
    % subplot(2,2,2)  → Altura del imán
    % subplot(2,2,3)  → Magnitud del campo magnético
    % subplot(2,2,4)  → Superficie 3D del campo

    %----------- ACTUALIZAR TIEMPO -----------%
    cc = cc + 1;                  % Avanza al siguiente paso de simulación
    t(cc) = t(cc-1) + dt;         % Avanza el tiempo
end
