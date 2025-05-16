%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Session 1 | Step 2    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rodrigo Gamboa & Francisco Montes | May 2024.


%After "Uncomment for setting a proper viewing perspective"
%view(-34,33) from previous step.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic-field Calculation & Visualization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Using the results obtained from the Biot-Savart law, get the B-field 
%at each point in space, around the current loop.


; ; ;   %Create three, 3x3 equal arrays, dBx,
                                           %dBy & dBz, each spanning from 1
                                           %to Lx, 1 to Ly and 1 to Lz, and 
                                           %initialize all values to zero. You should 
                                           %obtain in the workspace, three 101x101x101 
                                           % objects (dBx, dBy, dBz). 
                                           
                                           %What do you think we are to use
                                           %this for:________________________________?

                                    
for                       %Open a for loop using "I" as index, going from 1 to Lx.  
    for                   %Same but for Ly, use a "J" index.
        for               %Same but for Lz, use a "K" index.   

            for           %Open a 4th for loop, usin an "L" index, running from 1 to nl*N (Why?).

            %<inside the 4th loop>
                                
                                % Make a drawing (diagram), using a position vector 
                                % from the center of the coil to each current element, given by
                                % r´= Px(L)i+Py(L)j+Pz(L)k, and a position vector 
                                % from the center of the coil to the place in space  
                                % where we want to compute the B-fiel given by 
                                % R´= x(I)i+y(J)j+z(K)k. With this, find the position vector 
                                % from each current element to the place in space  
                                % where we want to compute the B-fiel, given by   
                                % r = (rx)i+(ry)j+(rz)k


            rx=;     %Write here the rx component
            ry=;     %Same but for the ry component
            rz=;     %Same but for the ry component

            r=;    %Get the magnitude

            ;                         %Declate an r3 variable
                                            %equal to the third power of r
                                           
            
            % Using Biot-Savart´s original expression, use a generic
            % dl = (dx)i+(dy)j and r = (rx)i+(ry)j+(rz)k to compute the 
            % cross product between dl and r, to find three independent 
            % expression for the components of the B-fiel we want to
            % compute, dBx(I,J,K), dBy(I,J,K) and dBz(I,J,K). Do not forget
            % to update or rewrite the former components, between each
            % cycle.

            ;              %Write here the dBx(I,J,K) component
            ;              %Same but for the dBy(I,J,K) component  
            ;   %Same but for the dBz(I,J,K) component

            %<inside the 4th loop>

            end   %Close the 4th loop

        end    %Close the 3rd loop   
    end   %Close the 2nd loop
end    %Close the 1st loop


%------------------------------------------------------------
%Uncomment all the following, and explain (comment) each line
%------------------------------------------------------------

%Bmag=sqrt(dBx.^2+dBy.^2+dBz.^2);    %?

%centery=round(Ly/2);                %?
%Bx_xz=squeeze(dBx(:,centery,:));    %?
%Bz_xz=squeeze(dBz(:,centery,:));    %?
%Bxz=squeeze(Bmag(:,centery,:));     %?

%figure(2)                           %?
%hold on

%pcolor(x,z,(Bxz').^(1/3)); shading interp; colormap jet; colorbar    %?
%h1=streamslice(x,z,Bx_xz',Bz_xz',3);                                 %?
%set(h1,'Color', [0.8 1 0.9]);                                        %?
%xlabel 'x'                                                           %? 
%ylabel 'z'                                                           %?
%title 'Magnetic field of a circular current'                         %?
%zz(length(ang))=0;                                                   %?
