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


%[Reserved word for definig a function] <function call>


%mo=4*pi*1e-7;           %Just uncomment.
%ds=0.005;               %Uncomment and continue ahead and then go back here to comment
                         %on the meaning of this variable.

                         ;      %Declare a x array, extending from -Rring to Rring in steps of ds.   
;      %Declare a y array, extending from -Rring to Rring in steps of ds.
;           %Define an array called Lx such that it has the
                         %same length than the array x.

;           %Define an array called Ly such that it has the
                         %same length than the array y.

;        %Define a 2D matrix, or array called Bz, such that
                         %it has entries from 1 to Lx and 1 to Ly, for the ith 
                         %and jth components, respectively, and initialize all 
                         %of them to zero.


;                 %Introduce a variable called phiB, and set it
                         %equal to zero.


%for i=1:Lx          %Start a for  loop using an "i" counter going from 1 to Lx.
    %for j=1:Ly            %Start a for  loop using a "j" counter going from 1 to Ly.

        %r=sqrt(x(i)^2+y(j)^2);       %Uncomment and comment what exactly is this?
                                      %MAKE A DRAWING!!!

        %if r<Rring                   %Why this? Comment and think!

            %Bz(i,j) = WRITE HERE THE EQUATION FOR THE FIELD FROM CANVAS
            % "STAGE 3: Calculation of induced fem by a magnet (dipole) free falling through a coil"

            %phiB=;   %Following intructions in CANVAS, write here the
                      %value for the B-Field Flux!

         %Close (end) the "if".   

    %Close the for with the "j" counter.  
%Close the for with the "i" counter.