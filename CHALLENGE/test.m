% Magnetic-field Calculation & Visualization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Array 3x3
dBx = zeros(Lx, Ly, Lz);% Componente X del campo magnético,
% inicializado en ceros
dBy = zeros(Lx, Ly, Lz);% Mismo que el anterior
% pero para y
dBz = zeros(Lx, Ly, Lz);% Mismo que el anterior
                        % pero para z

for I= 1:Lx %Un loop que va desde 1 a Lx 
   for J= 1:Ly %Mismo para y
       for K= 1:Lz %Igualmente pero para z
           for L= 1:nl*N %Otro for que va desde 1 a
               %el numero de loops multiplicado puntos por loop,
               %por si tuvieramos mas de un loop
            rp= Px(L)+Py(L)+Pz(L);%Sumamos las posiciones de x, y, z
            % en cada punto
            R= x(I)+y(J)+z(K);%Recorre todo el area de
            % trabajo y suma las coordenadas
           rx=x(I)-Px(L);%Componente r en x
           ry= y(J)-Py(L);%Mismo en y
           rz=z(K)-Pz(L);%Mismo en z
           r=sqrt(rx^2 + ry^2 + rz^2+ rw^2);%Esta es la magnitud de r
           r3= r^3;%Declaramos una variable r3 igual
           % a la magnitud de r elevado al cubo
           dBx(I,J,K) = dBx(I,J,K) + km * (dy(L)*rz) / r3;             
           %Definimos el diferencial del campo en x
           dBy(I,J,K) = dBy(I,J,K) - km * (dx(L)*rz) / r3;
           %Mismo que el pasado pero en y
           dBz(I,J,K) = dBz(I,J,K) + km * (dx(L)*ry - dy(L)*rx) / r3;
           %Mismo que los pasados pero en z
                                         
              
            end  
       end     
   end  
end

Bmag=sqrt(dBx.^2+dBy.^2+dBz.^2);    %Magnitud del
% Campo magnetico
centery=round(Ly/2);             % Round redondea
% la longitud de y entre 2
Bx_xz=squeeze(dBx(:,centery,:));    % Extrae la componente Bx del campo magnético en el plano XZ (fijando Y = centery)
Bz_xz = squeeze(dBz(:,centery,:));      % Extrae la componente Bz del campo magnético en el plano XZ (fijando Y = centery)
Bxz   = squeeze(Bmag(:,centery,:));     % Extrae la magnitud total del campo magnético en el plano XZ
figure(2)                           %Muestra la figura del campo
hold on
pcolor(x,z,(Bxz').^(1/3)); shading interp; colormap jet; colorbar   
% Muestra un mapa de color de la magnitud del campo magnético en el plano XZ
% Se eleva a la 1/3 para visualizar mejor valores bajos y altos
% 'shading interp' suaviza la interpolación entre colores
% 'colormap jet' aplica un degradado de colores tipo "arcoíris"
% 'colorbar' muestra la escala de colores a la derecha
h1=streamslice(x,z,Bx_xz',Bz_xz',3); % Dibuja líneas de
% flujo del campo usando las componentes Bx y Bz en el plano XZ
% El 3 controla la densidad de las líneas de flujo
%(entre más alto, menos densas)                       
set(h1, 'Color', [0.8 1 0.9]);          %Cambia el color de las líneas de flujo a un tono verde claro
xlabel('x')                              % Etiqueta del eje X
ylabel ('z')                             % Etiqueta del eje Z
title('Magnetic field of a circular current')  %Título de la gráfica
zz(length(ang))=0;



