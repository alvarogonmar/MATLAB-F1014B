%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Session 3 | Step 2              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rodrigo Gamboa & Francisco Montes | May 2024


%-----------------------------------------------------------------------
%Code continues from inside the while loop you used in the previous step
%-----------------------------------------------------------------------

 %<start  while>    %Start a while loop running until the position of the magnet (dinamically) in z 
                    %is bigger than -5 (i.e. the bottom of the canvas). 

                    %For diagnostics purposes, print the position of the magnet in z. 

    ; %Using the addpoints function mentioned above, add points (the zm(cc) ones) 
      %to the path, centered in x=0.

      %Use drawnow to update and/or modify graphics objects and 
      %want to see the updates on your canvas immediately.
    
    %head=scatter(0,zm(cc), 100,'filled'); %Uncomment this line and explain what its doing. 
   

    %%%%%%%%%%%%%%%%%
    %<NEW STUFF HERE>
    %%%%%%%%%%%%%%%%%

    %Go to the Canvas tab called "Stage 2: Magnetic Force between two magnetic
    %dipoles", read through it and perform the referred calculations. Using
    %this result. Get the force experienced by the Magnet Fm(cc), as a function
    %of z, i.e., as a function of it´s height. Use (adapt) your results to the 
    %variable names in yout code.


    ;                                        %Magnetic force Fm(cc)

    ;                                        %Net (total) force over the
                                             %magnet F(cc). Important, it 
                                             %is net force.
                                                                 
                                                                    

    ;      %Using the total force get the instantaneous acceleration and label it as "a".
           %Use the pause function, with a value of 0.01, to grasp the
           %head´s magnet movement.



                 ;    %Get new magnet´s position z(cc+1) using high school physics, 
                      %i.e. using constant __________________ kinematic equations. 

                 ;    %Do the same but for the free fall case,
                      %i.e., with no magnetic braking.

                 ;    %Get new magnet´s approx. velocity using basic velocity definition over
                      %a time interval dt, i.e., vz(cc+1) = delta(z)/delta(t).
    
                    %Just for fun, derive the exact analytic expression for the velocity from the expression 
                    %of the force. This is a Physics only step, however if you manage to incorporate it 
                    %into the code, and actually make it work, your team will get extra final points.    

                 ;  %Do the same (excepting the "just for fun steo") but for the free fall case,
                    %i.e., with no magnetic braking.
   

           ;                    %Increase the cc counter in one unit.

   %delete(head)                %Uncomment 
                                
   %%<end  while>             Close the while loop


%---------------------------------------------------------------------
%--Final Part---------------------------------------------------------
%---------------------------------------------------------------------

%Uncomment all of the following, and explain in 
%detail (in next step of the code, not now) what it does or will do. 

%figure(4)
%subplot(1,2,1)

%hold on
%plot(zm(1:length(Fm)),1000*Fm, '-b', 'LineWidth', 2) %Why do we have 1000Fm?
%plot([0,0],[-150,150],'-.k','LineWidth', 2)
%grid on
%xlabel 'z position (m)'
%ylabel 'Magnetic force (mN)'
%title 'Magnetic force of a Current ring over a falling magnet'
%legend('Magnetic force in the Z direction','Current loop location','Location','southeast')

%subplot(1,2,2)

%hold on
%tt=0:dt:(cc-1)*dt;
%plot(tt,zm,'-r', 'LineWidth', 2)
%plot(tt,zmfree,'--b', 'LineWidth', 2)   
%plot([0,1.8],[0,0],'-.k','LineWidth', 2)
%grid on
%xlabel 'time (s)'
%ylabel 'z position (m)'
%title 'Position vs time of a Magnetic dipole falling throug a current ring'
%legend('Fall over a current ring','Free fall (no Magnetic force)', 'Current loop location','Location','southwest')
%axis([0 1.8 -6 6])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINALLY: Run and test your code for limit cases!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%<THIS IS PART OF YOUR DELIVERABLE NO.1>


%(1) What happens if the magnetic moment is zero? (run and explain).

%(2) What happens if the mass of the magnet is set to 1 kg? (run and explain).

%(3) What happens if you increase the actual value for the magnetic
% moment one order if magnitude? (run and explain).

%(4) What happens if you decrease the mass of the magnet by one 
% order of magnitude? (run and explain).

%(5) How many plots are you getting at the end? Explain each in full
%detail. Full details please!



