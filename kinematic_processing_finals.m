function [Z_E, Z_F, Z_H, Z_P, Z_I] = kinematic_processing_finals(X,p)
%  this function 

%% UNPACKING INPUTS
    %geometry & linkage lengths
    el_EF = p.el_EF;
    el_HF = p.el_HF;
    el_FP = p.el_FP;
    el_G3P = p.el_G3P;
    el_G3H = (el_HF + el_FP) - el_G3P;
    
    %unpacking state vector = [th_HI; r_PI; thdot_HI; rdot_PI]
  % th_ED = X(:,1);
  % thdot_ED = X(:,2); 
    x_G1 = X(:,3);
    xdot_G1 = X(:,4);
    y_G1 = X(:,5);
    ydot_G1 = X(:,6);
    th_EF = X(:,7); 
    thdot_EF = X(:,8);
    x_G2 = X(:,9);
    xdot_G2 = X(:,10);
    y_G2 = X(:,11);
    ydot_G2 = X(:,12);
    th_HI = X(:,13);
    thdot_HI = X(:,14);
    x_G3 = X(:,15);
    xdot_G3 = X(:,16);
    y_G3 = X(:,17);
    ydot_G3 = X(:,18);
    r_PI = X(:,19);
    rdot_PI = X(:,20);
    
%% KINEMATIC EQUATIONS (JOINT POSITIONS & VELOCITIES)
    %linkage point E kinematics
    x_E = 2.*x_G1;                  
    y_E = 2.*y_G1;              
    xdot_E = 2.*xdot_G1;
    ydot_E = 2.*ydot_G1;
    
    %linkage point F kinematics
    x_F = x_G2 - el_EF.*sin(th_EF)./2;
    y_F = y_G2 + el_EF.*cos(th_EF)./2;
    xdot_F = xdot_G2 - thdot_EF.*el_EF.*cos(th_EF)/2;
    ydot_F = ydot_G2 - thdot_EF.*el_EF.*sin(th_EF)/2;
    
    %linkage point P kinematics
    x_P = x_G3 + el_G3P.*sin(th_HI); 
    y_P = y_G3 - el_G3P.*cos(th_HI);
    xdot_P = xdot_G3 + thdot_HI.*el_G3P.*cos(th_HI);
    ydot_P = ydot_G3 + thdot_HI.*el_G3P.*sin(th_HI);
    
    %linkage point H kinematics
    x_H = x_G3 - el_G3H.*sin(th_HI); 
    y_H = y_G3 + el_G3H.*cos(th_HI);
    xdot_H = xdot_G3 - thdot_HI.*el_G3H.*cos(th_HI);
    ydot_H = ydot_G3 - thdot_HI.*el_G3H.*sin(th_HI);
    
    %linkage point I kinematics
    x_I = x_P + r_PI.*sin(th_HI);
    y_I = y_P - r_PI.*cos(th_HI);
    xdot_I = xdot_P + rdot_PI.*sin(th_HI) + r_PI.*thdot_HI.*cos(th_HI);
    ydot_I = ydot_P - rdot_PI.*cos(th_HI) + r_PI.*thdot_HI.*sin(th_HI);
    
%% OUTPUT VECTORS
    Z_E = [x_E, y_E, xdot_E, ydot_E];   %position & velocity of point E
    Z_F = [x_F, y_F, xdot_F, ydot_F];   %position & velocity of point F
    Z_H = [x_H, y_H, xdot_H, ydot_H];   %position & velocity of point H
    Z_P = [x_P, y_P, xdot_P, ydot_P];   %position & velocity of point P
    Z_I = [x_I, y_I, xdot_I, ydot_I];   %position & velocity of point I
     
end
