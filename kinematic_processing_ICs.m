function [A_G, A_joint, Z_ang, Z_leng] = kinematic_processing_ICs(Y, thdot_ED, R, p)
%% FUNCTION DESCRIPTION
% This function uses kinematic constraints on the linkage system to get
% four main groups of outputs from just 3 input angles (1DOF). It
% first determines all the positions by altering point P's, then uses a
% matrix equation to simultaneously solve all the velocity equations
% (since we cannot alter the velocity of point P being 0) and gives
% the resulting angular velocities for each linkage. 
% 
% This function therefore allows the user to define any set of angles 
% as their initial condition that they want, while still having a valid 
% dynamic system, despite not every initial condition BEING valid for a 
% PREDEFINED set of geometry. 
% 
% This function can ALSO be used in a practical sense by 
% iterating through every possible Y input in order to find what 
% P-position you would like for the robot to start in.
% 
%% INPUT/OUTPUT DESCRIPTION
% INPUTS:
%  1: The desired initial angle state Y, which specifies the desired angle 
%     (CCW) of each linkage with respect to the +jhat direction
%      where Y = [th_ED; th_EF; th_HI]
%  2: the desired input angle velocity ED (we can still specify this
%     because it is a part of our 1 DOF that we still have for velocity
%  3: input parameters, these are just the predefined constants for the
%     problem
%  4: end effector parameters
%     where R = [r_PI; rdot_PI]
% 
% OUTPUTS:
%  1: the position and velocity of each linkage's C.O.M.
%     (Z_G1, Z_G2, Z_G3) where Z = [x; xdot; y; ydot]
%
%  2: the position and velocity of each pin joint in the system
%     (Z_E, Z_F, Z_H, Z_P, Z_I) where Z = [x; xdot; y; ydot]
%
%  3: the angular velocity (CCW) of each linkage with respect to the 
%     +jhat direction
%     (Z_ang) where Z = [thdot_ED; thdot_EF; thdot_HI]
%  
%  4: the vertical and horizontal position of point P, as well as the
%     lengths of 
%     this allows for TWO extra specified angles, theta_HI and theta_EF, to
%     be input as a desired value, and this function will then give the
%     two lengths h1 and d1. This effectively allows us to take the
%     1 DOF system, and add an extra 2 DOF FOR THE INPUT only, via changing
%     the geometry of the problem statement. This method is motivated 
%     by the nonlinear dynamics which would have been required to solve
%     otherwise (using a predefined P position).

%% UNPACKING INPUTS
    %predefined linkage length parameters
    el_ED = p.el_ED;
    el_EF = p.el_EF;
    el_HF = p.el_HF;
    el_FP = p.el_FP;
    el_G3F = (el_HF - el_FP)/2;
    
    %input vector
    th_ED = Y(1);
    th_EF = Y(2);
    th_HI = Y(3);
    
    r_PI = R(1);
    rdot_PI = R(2);
    
%% KINEMATIC EQUATIONS (JOINT POSITIONS)
    %linkage point E kinematics
    x_E = -el_ED*sin(th_ED);
    y_E = el_ED*cos(th_ED);

    %linkage point F kinematics
    x_F = x_E - el_EF*sin(th_EF);
    y_F = y_E + el_EF*cos(th_EF);
    
    %linkage point H kinematics
    x_H = x_F - el_HF*sin(th_HI);
    y_H = y_F + el_HF*cos(th_HI);
    
    %determining the position of point P
    x_P = x_F + el_FP*sin(th_HI); 
    y_P = y_F - el_FP*cos(th_HI);
    d1 = x_P;
    h1 = -y_P;
    
    %linkage point I kinematics
    x_I = x_P + r_PI*sin(th_HI);
    y_I = y_P - r_PI*cos(th_HI);
    
%% KINEMATIC EQUATIONS (C.O.M. POSITIONS)
    %center of mass G1 kinematics
    x_G1 = -el_ED*sin(th_ED)/2;
    y_G1 = el_ED*cos(th_ED)/2;

    %center of mass G2 kinematics
    x_G2 = x_E - el_EF*sin(th_EF)/2;
    y_G2 = y_E + el_EF*cos(th_EF)/2;
    
    %center of mass G3 kinematics
    x_G3 = x_F - el_G3F*sin(th_HI); 
    y_G3 = y_F + el_G3F*cos(th_HI);

%% KINEMATIC EQUATIONS (JOINT VELOCITIES)
    A_vel = [1, 0, el_ED*cos(th_ED), zeros(1,6);
             0, 1, el_ED*sin(th_ED), zeros(1,6);
             -1, 0, 0, 1, 0, el_EF*cos(th_EF), 0, 0, 0;
             0, -1, 0, 0, 1, el_EF*sin(th_EF), 0, 0, 0;
             0, 0, 0, -1, 0, 0, 1, 0, -el_FP*cos(th_HI);
             0, 0, 0, 0, -1, 0, 0, 1, -el_FP*sin(th_HI);
             zeros(1,6), 1, 0, 0;
             zeros(1,7), 1, 0;
             0, 0, 1, zeros(1,6)];
             
    B_vel = [zeros(8,1);
             thdot_ED];
         
    W = A_vel\B_vel;
    
    xdot_E = W(1);
    ydot_E = W(2);
  % thdot_ED = W(3); %already defined in the inputs
    xdot_F = W(4);
    ydot_F = W(5);
    thdot_EF = W(6);
    xdot_P = W(7);
    ydot_P = W(8);
    thdot_HI = W(9);
    
    %extra kinematic equations which are just extensions of point P
    xdot_H = xdot_F - thdot_HI*el_HF*cos(th_HI);
    ydot_H = ydot_F - thdot_HI*el_HF*sin(th_HI);
    xdot_I = xdot_P + rdot_PI*sin(th_HI) + r_PI*thdot_HI*cos(th_HI);
    ydot_I = ydot_P - rdot_PI*cos(th_HI) + r_PI*thdot_HI*sin(th_HI);

%% KINEMATIC EQUATIONS (C.O.M. VELOCITIES)
    %center of mass G1 kinematics
    xdot_G1 = -thdot_ED*el_ED*cos(th_ED)/2;
    ydot_G1 = -thdot_ED*el_ED*sin(th_ED)/2;
    
    %center of mass G2 kinematics
    xdot_G2 = xdot_E - thdot_EF*el_EF*cos(th_EF)/2;
    ydot_G2 = ydot_E - thdot_EF*el_EF*sin(th_EF)/2;
    
    %center of mass G3 kinematics
    xdot_G3 = xdot_F - el_G3F*thdot_HI*cos(th_HI); 
    ydot_G3 = ydot_F - el_G3F*thdot_HI*sin(th_HI);
    
%% OUTPUT VECTORS 
    %C.O.M. outputs
    Z_G1 = [x_G1;       %position & velocity of the C.O.M. for linkage ED
            y_G1;
            xdot_G1;
            ydot_G1];
    Z_G2 = [x_G2;       %position & velocity of the C.O.M. for linkage EF
            y_G2;
            xdot_G2;
            ydot_G2];
    Z_G3 = [x_G3;       %position & velocity of the C.O.M. for linkage HP
            y_G3;
            xdot_G3;
            ydot_G3];
    A_G = [Z_G1, Z_G2, Z_G3]; %packaging vectors into a matrix
       
    %pin joint outputs
    Z_E = [x_E;         %position & velocity of point E
           y_E;
           xdot_E;
           ydot_E];
    Z_F = [x_F;         %position & velocity of point F
           y_F;
           xdot_F;
           ydot_F]; 
    Z_H = [x_H;         %position & velocity of point H
           y_H;
           xdot_H;
           ydot_H];
    Z_P = [x_P;         %position & velocity of point P
           y_P;
           xdot_P;
           ydot_P];
    Z_I = [x_I;         %position & velocity of point I
           y_I;
           xdot_I;
           ydot_I];
       
    A_joint = [Z_E, Z_F, Z_H, Z_P, Z_I]; %packaging vectors into a matrix
       
    %angle output
    Z_ang = [thdot_ED;  %linkage angles & angular velocities
             thdot_EF;
             thdot_HI];

    %geometry lengths output
    Z_leng = [h1;       %horizontal and vertical position of point P
              d1];
end
