mag= 2000;        %(Magnetic moment of the magnet)
                 %Intensidad y orientacion del "Dipolo Magentico
                 %and initialize to 2000. What are its units? "A*m^2"


M_mass= 0.004;    %Declare the Magnet´s mass in kg, use 0.004 to begin.
                 %How many grams is this= 0.004kg---> 4 gr?
g= 9.81;          %Gravedad
w=M_mass*g;       %Magnet´s weigth in N ---> Masa x gravedad = peso

zo=5;             %Magnet's initial position


dt=0.05;          %time step 


zm(1)=zo;         %(magnet position) and store in its first value the Magnet's initial position
                 %Cuando hay Fuerzas MAGNÉTICAS actuando


zmfree(1)=zo;    %"zmfree", which stands for the magnte free fall case
                %think we introduced this --> Guardar las posiciones de caida libre
                %cuando no hay fuerzas actuando aparte de la gravedad


tt(1)=0;         %Declare a vector called "tt" (time) vector de tiempo


vz(1)=0;         %Declare a vector called "vz" (Z component for velocity)


vzfree(1)=0;     % "vzfree". Again, why? -> Obtener la velocidades de mi magneto en caso de caida libre sin fuerzas
                % aparte de la gravedad


cc=1;            %cc=1 --> Contador del ciclo


path=animatedline;      %Use the animatedline function (use Matlab´s upper-right search function bar) to create an
                       %animated line that has no data and adds it to the current axes. Later (a few lines below),
                       %add points to the line in a loop using the addpoints function, add data.


%<start  while>        %Start a while loop running until the position of the magnet (dinamically) in z
while zm(cc) > -5
     zm(cc)              %For diagnostics purposes, print the position of the magnet in z. 

    addpoints(h, 0, zm(cc));   %Using the addpoints function mentioned above, add points (the zm(cc) ones) 
                    %to the path, centered in x=0.
    
    
    drawnow               %Use drawnow to update and/or modify graphics objects and 
                    %want to see the updates on your canvas immediately.
    
head=scatter(0,zm(cc), 100,'filled');       %Uncomment this line and explain what its doing. 
   
    %%%%%%%%%%%%%%%%%
    %<NEW STUFF HERE>
    %%%%%%%%%%%%%%%%%

    %Go to the Canvas tab called "Stage 2: Magnetic Force between two magnetic
    %dipoles", read through it and perform the referred calculations. Using
    %this result. Get the force experienced by the Magnet Fm(cc), as a function
    %of z, i.e., as a function of it´s height. Use (adapt) your results to the 
    %variable names in yout code.


    Fm(cc) = -mag * (mo * I * R^2 * (-3) * zm(cc)) / (2 * (zm(cc)^2 + R^2)^(5/2));             %Magnetic force Fm(cc)

    F(cc) = -w + Fm(cc);                     %Net (total) force over the
                                             %magnet F(cc). Important, it 
                                             %is net force.
                                                                 
                                                                    

    a = F(cc) / mag_mass;     %Using the total force get the instantaneous acceleration and label it as "a".
    pause(0.01);       %Use the pause function, with a value of 0.01, to grasp the
                        %head´s magnet movement.



    zm(cc+1) = zm(cc) + vz(cc)*dt + 0.5*a*dt^2;        %Get new magnet´s position z(cc+1) using high school physics, 
                                                           %i.e. using constant __________________ kinematic equations. 

    zmfree(cc+1) = zmfree(cc) + vzfree(cc)*dt - 0.5*9.81*dt^2;               %Do the same but for the free fall case,
                      %i.e., with no magnetic braking.

    vz(cc+1) = (zm(cc+1) - zm(cc)) / dt;         %Get new magnet´s approx. velocity using basic velocity definition over
                                                  %a time interval dt, i.e., vz(cc+1) = delta(z)/delta(t).
    
                    %Just for fun, derive the exact analytic expression for the velocity from the expression 
                    %of the force. This is a Physics only step, however if you manage to incorporate it 
                    %into the code, and actually make it work, your team will get extra final points.    

    vzfree(cc+1) = (zmfree(cc+1) - zmfree(cc)) / dt;           %Do the same (excepting the "just for fun steo") but for the free fall case,
                                                                %i.e., with no magnetic braking.
   

     cc = cc + 1;                    %Increase the cc counter in one unit.

   delete(head)                %Uncomment 
                                
   %%<end  while>             Close the while loop

    
    %delete(head)                %Leave this commented.


end %             Close the while loop

%---------------------------------------------------------------------
%--Final Part---------------------------------------------------------
%---------------------------------------------------------------------

%Uncomment all of the following, and explain in 
%detail (in next step of the code, not now) what it does or will do. 

figure(4)
subplot(1,2,1)

hold on
plot(zm(1:length(Fm)),1000*Fm, '-b', 'LineWidth', 2) %Why do we have 1000Fm? -> Transaforma
                                                    % Newtons a miliNewtons 1N=1000mN
plot([0,0],[-150,150],'-.k','LineWidth', 2)
grid on
xlabel 'z position (m)'
ylabel 'Magnetic force (mN)'
title 'Magnetic force of a Current ring over a falling magnet'
legend('Magnetic force in the Z direction','Current loop location','Location','southeast')

subplot(1,2,2)

hold on
tt=0:dt:(cc-1)*dt;
plot(tt,zm,'-r', 'LineWidth', 2)
plot(tt,zmfree,'--b', 'LineWidth', 2)   
plot([0,1.8],[0,0],'-.k','LineWidth', 2)
grid on
xlabel 'time (s)'
ylabel 'z position (m)'
title 'Position vs time of a Magnetic dipole falling throug a current ring'
legend('Fall over a current ring','Free fall (no Magnetic force)', 'Current loop location','Location','southwest')
axis([0 1.8 -6 6])

%---------------------------------------------------------------------
%--Final Part---------------------------------------------------------
%---------------------------------------------------------------------

%Uncomment all of the following, and explain in 
%detail (in next step of the code, not now) what it does or will do. 

figure(4)
subplot(1,2,1)

hold on
plot(zm(1:length(Fm)),1000*Fm, '-b', 'LineWidth', 2) %Why do we have 1000Fm?
plot([0,0],[-150,150],'-.k','LineWidth', 2)
grid on
xlabel 'z position (m)'
ylabel 'Magnetic force (mN)'
title 'Magnetic force of a Current ring over a falling magnet'
legend('Magnetic force in the Z direction','Current loop location','Location','southeast')

subplot(1,2,2)

hold on
tt=0:dt:(cc-1)*dt;
plot(tt,zm,'-r', 'LineWidth', 2)
plot(tt,zmfree,'--b', 'LineWidth', 2)   
plot([0,1.8],[0,0],'-.k','LineWidth', 2)
grid on
xlabel 'time (s)'
ylabel 'z position (m)'
title 'Position vs time of a Magnetic dipole falling throug a current ring'
legend('Fall over a current ring','Free fall (no Magnetic force)', 'Current loop location','Location','southwest')
axis([0 1.8 -6 6])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINALLY: Run and test your code for limit cases!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<THIS IS PART OF YOUR DELIVERABLE NO.1>


%(1) What happens if the magnetic moment is zero? (run and explain).
% if the magnetic moment is zero the magnet falss only with the gravity
% force, just like a normal object in free fall.

%(2) What happens if the mass of the magnet is set to 1 kg? (run and explain).
% If the magnet’s mass is increased to 1 kg, the magnet will fall much
% faster,  because the heavier the magnet, the less effect the magnetic force has on stopping it.

%(3) What happens if you increase the actual value for the magnetic
% moment one order if magnitude? (run and explain).
% The magnet will be much more affected by the coil, so it will fall more
% slowly or depending wich value you define, is not going to fall, its
% going to be bouncing

%(4) What happens if you decrease the mass of the magnet by one 
% order of magnitude? (run and explain).
% The magnet becomes much lighter, the magnet will fall much slower, it
% takes a lot longer to reach the bottom, and again its going to be
% bouncing

%(5) How many plots are you getting at the end? Explain each in full
%detail. Full details please!
% We have two plots in the leff we have "Magnetic force of a Current ring over a falling magnet"
% this plot shows how the magnetic force acts on the falling magnet as it moves up and down along the z-axis.
% The blue curve shows how the force changes as the magnet falls. 
% and when the magnet is near the coil (z = 0), the magnetic force is strongest.
% and the force changes direction depending on where the magnet is.

% And in the right plot we have the "Position vs time of a Magnetic dipole falling through a current ring"
% this plot shows how the magnet’s position changes over time as it falls,
% in the horizonta axis we have time, and in the vertical axis we have z
% position, we have a red curve that is the path of the magnet falling with
% the effect of the magnetic force, and we have blue curve that is the path of the magnet 
% in free fall, with no magnetic force, just using gravity.