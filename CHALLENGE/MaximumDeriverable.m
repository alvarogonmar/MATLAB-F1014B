%---------------------------------------------------------%
%     Simulación de imán cayendo por una espira con FEM   %
%---------------------------------------------------------%
clc; clear; close all;

%------------------------ ESPIRA ------------------------%
nl = 1;
N = 20;
totalPoints = N * nl;

ds = 0.1;
x = -5:ds:5;
y = -5:ds:5;
z = x;
Lx = length(x);
Ly = length(y);
Lz = length(z);

rw = 0.2;
I = 300;
mo = 4*pi*1e-7;
km = mo*I/(4*pi);
R = 1.5;
dtheta = 2*pi/N;
ang = 0:dtheta:(2*pi - dtheta);

Px = zeros(1, totalPoints);
Py = zeros(1, totalPoints);
Pz = zeros(1, totalPoints);
dx = zeros(1, totalPoints);
dy = zeros(1, totalPoints);
dz = zeros(1, totalPoints);

s = 1;
for I = 1:nl
    Px(s:s+N-1) = R * cos(ang);
    Py(s:s+N-1) = R * sin(ang);
    Pz(s:s+N-1) = 0;  % Coloca la espira en z=0
    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;
    dy(s:s+N-1) = Px(s:s+N-1) * dtheta;
    s = s + N;
end

dz = zeros(1, totalPoints);

%-------------------- PARÁMETROS DEL IMÁN --------------------%
mag = 500;
Rring = 0.5;
zo = 0.1;
dt = 0.01;
zring = 0;

t = zeros(1, 2000);
zm = zeros(1, 2000);
zm(1) = zo;
vz = zeros(1, 2000);
fem = zeros(1, 2000);
cc = 1;

%------------------ SIMULACIÓN PRINCIPAL ------------------%
figure(1); % Para gráficas fem, posición y campo

while zm(cc) > 0.0172
   
    pause(0.001);

    %----------- ANIMACIÓN EN FIGURA 99 -----------%
    figure(99); clf
    quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2)
    hold on
    scatter3(0, 0, zm(cc), 100, 'filled')
    axis([-2 2 -2 2 -0.2 0.2])
    xlabel('x'); ylabel('y'); zlabel('z');
    title(['Posición del imán (z = ' num2str(zm(cc), '%.3f') ' m)'])
    view(3)
    drawnow

    %----------- CÁLCULO FEM -----------%
    [xg, yg, phiB1, Bz] = B_due_M(zm(cc), mag, Rring);
    zm(cc+1) = zm(cc) + vz(cc)*dt - 0.5*9.81*dt^2;
    vz(cc+1) = vz(cc) - 9.81*dt;
    [~, ~, phiB2, ~] = B_due_M(zm(cc+1), mag, Rring);
    fem(cc) = - (phiB2 - phiB1) / dt;

    %----------- GRAFICAR FEM Y CAMPO EN FIGURA 1 -----------%
    figure(1);
    subplot(2,2,1)
    hold on; grid on
    xlabel('time, s'); ylabel('fem, mV')
    plot(t(1:cc), 100*fem(1:cc), '-k')
    plot(t(1:cc), 100*fem(1:cc), '*r')

    subplot(2,2,2)
    hold on; grid on
    xlabel('time, s'); ylabel('magnet height, cm')
    plot(t(1:cc), 100*zm(1:cc), 'ob')

    subplot(2,2,3)
    hold on
    pcolor(xg, yg, zm(cc)/abs(zm(cc)) * abs(0.005^2 * Bz).^(1/3));
    shading interp; colormap hot; colorbar
    view(-45,-45)

    subplot(2,2,4)
    hold on
    mesh(xg, yg, zm(cc)/abs(zm(cc)) * abs(10^2 * Bz).^(1/3));
    view(-30, -3)
    axis([-Rring Rring -Rring Rring -5 15])

    % Actualizar contador
    cc = cc + 1;
    t(cc) = t(cc-1) + dt;
end