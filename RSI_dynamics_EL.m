function [Xdot] = RSI_dynamics_EL(t,X,p)
%  this function creates the EOM for the robotic surgery dynamic system
%  using the Euler-Lagrange (EL) approach

%unpacking parameter quantities
    %linkage moments of inertia
    I_ED = p.I_ED;  
    I_EF = p.I_EF;
    I_HP = p.I_HP;

    %geometry & linkage lengths
    h1 = p.h1;    
    d1 = p.d1;
    el_ED = p.el_ED;
    el_EF = p.el_EF;
    el_HF = p.el_HF;
    el_FP = p.el_FP;
    el_G3F = (el_HF - el_FP)/2;

    %linkage masses
    m_ED = p.m_ED;  
    m_EF = p.m_EF;
    m_HP = p.m_HP;

%kinematic processing
    [Z_E, Z_F, ~, Z_P, Z_I] = kinematic_processing(X,p);
    
    %point E
    x_E = Z_E(1);
    y_E = Z_E(2);
    xdot_E = Z_E(3);
    ydot_E = Z_E(4);
    
    %point F
    x_F = Z_F(1);
    y_F = Z_F(2);
    xdot_F = Z_F(3);
    ydot_F = Z_F(4);
    
    %point P
    x_P = Z_P(1);
    y_P = Z_P(2);
    xdot_P = Z_P(3);
    ydot_P = Z_P(4);
        
    %point I
    x_I = Z_I(1);
    y_I = Z_I(2);
    xdot_I = Z_I(3);
    ydot_I = Z_I(4);
    
%Euler-Lagrange equations:

%Output: derivative of the state vector
    Xdot = [X(3);
            X(4);
            thddot_HI;
            rddot_PI];
end

