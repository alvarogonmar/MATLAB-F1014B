%-------------------------------%
%   Challenge Session No. 4
%-------------------------------%
%Rodrigo Gamboa & Francisco Montes | May 2024


%--------Part 4---------% fem CALCULATION | DECLARATION
    
    %-- After the second instance call for 
    % [x,y,phiB2,Bz]=B_due_M(zm(cc+1),mag,Rring);
     
  
    %Taking into account that the B_due_M function has been called 
    %twice (why)? Calculate a vector function of fem(cc) equal to the
    %"rate change or difference of the magnetic flux, with respect to
    % time", associated with the area of the ring (coil or wire loop).

    %fem(cc)=?

    %Uncomment all of the following, and explain in detals what each line 
    %is doing.

    %subplot(2,2,1)
    %hold on
    %grid on
    %xlabel 'time, s'
    %ylabel 'fem, mV'
    %plot(t(1:cc),100*fem(1:cc),'-k','LineWidth',1)
    %plot(t(1:cc),100*fem(1:cc),'*r','LineWidth',2)
    

    %subplot(2,2,2)
    %hold on
    %axis([0 0.3 -10 10])
    %grid on
    %xlabel 'time, s'
    %ylabel 'magnet heigth, cm'
    %plot(t(1:cc),100*zm(1:cc),'ob','LineWidth',2)    

    %subplot(2,2,3)
    %hold on
    %pcolor(x,y,zm(cc)/abs(zm(cc))*abs(abs(0.005^2*Bz)).^(1/3)); shading interp; colormap hot; colorbar
    %view(-45,-45)

    %subplot(2,2,4)
    
    %hold on
    %mesh(x,y,zm(cc)/abs(zm(cc))*abs(abs(10^2*Bz)).^(1/3)); 
    %view(-30,-3)
    %axis([-Rring Rring -Rring Rring -5 15])
    
    %cc=cc+1;
    %t(cc)=t(cc-1)+dt;
   
% close the while loop

%END OF CODE !!!%




