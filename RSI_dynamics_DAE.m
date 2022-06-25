function [Xdot] = RSI_dynamics_DAE(t,X,p)
%% FUNCTION DESCRIPTION []
%  this function creates the EOM for the robotic surgery dynamic system
%  using the Differential Algebraic Equations (DAE) approach

%% UNPACKING PREDEFINED PARAMETERS []
    %gravitational acceleration
    g = p.g; 

    %linkage moments of inertia
    I_ED = p.I_ED;  
    I_EF = p.I_EF;
    I_HP = p.I_HP;

    %geometry & linkage lengths
    el_ED = p.el_ED;
    el_EF = p.el_EF;
    el_HF = p.el_HF;
    el_FP = p.el_FP;
%     el_G3P = p.el_G3P;
    el_G3F = (el_HF - el_FP)/2;

    %linkage masses
    m_ED = p.m_ED;  
    m_EF = p.m_EF;
    m_HP = p.m_HP;

%% STATE VECTOR PROCESSING []
    th_ED = X(1);       %link ED state
    thdot_ED = X(2); 
    x_G1 = X(3); 
    xdot_G1 = X(4);
    y_G1 = X(5);
    ydot_G1 = X(6);
    th_EF = X(7);       %link EF state
    thdot_EF = X(8); 
  % x_G2 = X(9); 
    xdot_G2 = X(10); 
  % y_G2 = X(11); 
    ydot_G2 = X(12);
    th_HI = X(13);      %link HI state
    thdot_HI = X(14); 
  % x_G3 = X(15); 
    xdot_G3 = X(16); 
  % y_G3 = X(17); 
    ydot_G3 = X(18);
  % r_PI = X(19);       %end effector state
    rdot_PI = X(20);
  % x_P = X(21);        %point P state
    xdot_P = X(22);
  % y_P = X(23);
    ydot_P = X(24);
  % x_F = X(25);        %point F state
    xdot_F = X(26);
  % y_F = X(27);
    ydot_F = X(28);
   
%% KINEMATIC PROCESSING
    %miscellaneous parameters
    S = el_EF*sin(th_EF)/2;
    C = el_EF*cos(th_EF)/2;
    S2 = (el_HF - el_FP)*sin(th_HI)/2;
    C2 = (el_HF - el_FP)*cos(th_HI)/2;
    
    %motor equations
    Mmot_D = 0;
    Mmot_E = 0;
    Mmot_F = 0;
    
    %end effector
    omega = 2*pi;            %freq that end effector oscillates at [rad/s]
    rddot_PI = omega*cos(omega*t);
        
%% DAE MATRICES
 A = [y_G1, -x_G1, y_G1, -x_G1, 0, 0, -I_ED, zeros(1,12);       %eq(1)
     1, 0, -1, 0, 0, 0, 0, -m_ED, zeros(1,11);                  %eq(2)
     0, 1, 0, -1, 0, 0, 0, 0, -m_ED, zeros(1,10);               %eq(3)
     0, 0, 1, 0, -1, zeros(1,5), -m_EF, zeros(1,8);             %eq(4)
     0, 0, 0, 1, 0, -1, zeros(1,5), -m_EF, zeros(1,7);          %eq(5)
     0, 0, C, S, C, S, 0, 0, 0, -I_EF, zeros(1,9);              %eq(6)
     0, 0, 0, 0, 1, zeros(1,8), -m_HP, zeros(1,5);              %eq(7)
     zeros(1,5), 1, zeros(1,8), -m_HP, 0, 0, 0, 0;              %eq(8)
     0, 0, 0, 0, C2, S2, zeros(1,6), -I_HP, zeros(1,6);         %eq(9)
     zeros(1,6), el_ED*cos(th_ED)/2, 1, zeros(1,11);            %eq(b1)
     zeros(1,6), el_ED*sin(th_ED)/2, 0, 1, zeros(1,10);         %eq(b2)
     zeros(1,7), -2, 0, el_EF*cos(th_EF), zeros(1,7), 1, 0;      %eq(b3)
     zeros(1,8), -2, el_EF*sin(th_EF), zeros(1,8), 1;      %eq(b4)
     zeros(1,12), el_G3F*cos(th_HI), 1, 0, 0, 0, -1, 0;         %eq(b5)
     zeros(1,12), el_G3F*sin(th_HI), 0, 1, 0, 0, 0, -1;         %eq(b6)
     zeros(1,15), 1, 0, 0, 0;                                   %eq(b7)
     zeros(1,15), 0, 1, 0, 0;                                   %eq(b8)
     zeros(1,12), -el_FP*cos(th_HI), 0, 0, 1, 0, -1, 0;         %eq(b9)
     zeros(1,12), -el_FP*sin(th_HI), 0, 0, 0, 1, 0, -1];        %eq(b10)
 
B = [Mmot_E - Mmot_D;                                           %eq(1)
     0;                                                         %eq(2)
     m_ED*g;                                                    %eq(3)
     0;                                                         %eq(4)
     m_EF*g;                                                    %eq(5)
     Mmot_F - Mmot_E;                                           %eq(6)
     0;                                                         %eq(7)
     m_HP*g;                                                    %eq(8)
     -Mmot_F;                                                   %eq(9)
     el_ED*thdot_ED^2*sin(th_ED)/2;                             %eq(b1)
     -el_ED*thdot_ED^2*cos(th_ED)/2;                            %eq(b2)
     el_EF*thdot_EF^2*sin(th_EF);                               %eq(b3)
     -el_EF*thdot_EF^2*cos(th_EF);                              %eq(b4)
     el_G3F*thdot_HI^2*sin(th_HI);                              %eq(b5)
     -el_G3F*thdot_HI^2*cos(th_HI);                             %eq(b6)
     0;                                                         %eq(b7)
     0;                                                         %eq(b8)
     -el_FP*thdot_HI^2*sin(th_HI);                              %eq(b9)
     el_FP*thdot_HI^2*cos(th_HI)];                              %eq(b10)
     
%% SOLVING SYSTEM OF LINEAR EQUATIONS
    Z = A\B;
    thddot_ED = Z(7);
    xddot_G1 = Z(8);
    yddot_G1 = Z(9);
    thddot_EF = Z(10);
    xddot_G2 = Z(11);
    yddot_G2 = Z(12);
    thddot_HI = Z(13);
    xddot_G3 = Z(14);
    yddot_G3 = Z(15);
    xddot_P = Z(16);
    yddot_P = Z(17);
    xddot_F = Z(18);
    yddot_F = Z(19);
    
%% OUTPUT (DERIVATIVE OF STATE VECTOR)
    Xdot = [thdot_ED;   %link ED state
            thddot_ED;
            xdot_G1;
            xddot_G1;
            ydot_G1;
            yddot_G1;   
            thdot_EF;   %link EF state
            thddot_EF;
            xdot_G2;
            xddot_G2;
            ydot_G2;
            yddot_G2;
            thdot_HI;   %link HI state
            thddot_HI;
            xdot_G3;
            xddot_G3;
            ydot_G3;
            yddot_G3;
            rdot_PI;    %end effector state
            rddot_PI
            xdot_P;     %point P state
            xddot_P;
            ydot_P;
            yddot_P
            xdot_F;     %point F state
            xddot_F
            ydot_F
            yddot_F];
end