%% Final Project - Robotic Surgery Instrument []
% Created by Sean Hughes
% Date: 12/09/21
    clc; clear; close all;
    
    EL = false;
    DAE = true;
    
%% Givens / Constants []
    p = [];

    %gravitational acceleration [m/s^2]
    p.g = 9.81; 

    %geometry & linkage lengths [m]
    p.el_ED = 13*(0.0254);
    p.el_EF = 15*(0.0254);
    p.el_HF = 25*(0.0254);
    p.el_FP = 10*(0.0254);
    el_HP = p.el_HF + p.el_FP;
    p.el_G3P = el_HP/2;

    %linkage masses [kg]
    p.m_ED = 50*(0.453592);  
    p.m_EF = 55*(0.453592);
    p.m_HP = 80*(0.453592);
    m_PI = 0;

    %linkage moments of inertia [kg*m^2]
    p.I_ED = p.m_ED*p.el_ED^2/12;  
    p.I_EF = p.m_EF*p.el_EF^2/12;
    p.I_HP = p.m_HP*el_HP^2/12;

%% Initial Conditions []
    %state vector for EL 
    X0_EL = [  pi/4;  1;  0;  0];   %X0 = [th_HI; r_PI; thdot_HI; rdot_PI]
%     X0_EL = [ -pi/4;  1;  0;  0];
%     X0_EL = [     0;  5;  1; -1];
%     X0_EL = [3*pi/8;  0; -1;  1];
    
    % state vector for DAE
    Y0 = [0; 3*pi/2; pi/4];    %Y0 = [th_ED; th_EF; th_HI]
    thdot_ED0 = 0;
    R0 = [0.1; 0];
    
    [A_G, A_joint, Z_ang, Z_leng] = kinematic_processing_ICs(Y0, thdot_ED0, R0, p);
    
    %X0     = [th_ED; thdot_ED; x_G1; xdot_G1; y_G1; ydot_G1;
    %          th_EF; thdot_EF; x_G2; xdot_G2; y_G2; ydot_G2;
    %          th_HI; thdot_HI; x_G3; xdot_G3; y_G3; ydot_G3;
    %          r_PI; rdot_PI; 
    %          x_P; xdot_P; y_P; y_dot_P;
    %          x_F; xdot_F; y_F; y_dot_F];
    X0_DAE = [Y0(1); Z_ang(1); A_G(1,1); A_G(3,1); A_G(2,1); A_G(4,1);
              Y0(2); Z_ang(2); A_G(1,2); A_G(3,2); A_G(2,2); A_G(4,2);
              Y0(3); Z_ang(3); A_G(1,3); A_G(3,3); A_G(2,3); A_G(4,3);
              R0(1); R0(2);
              A_joint(1,4); A_joint(3,4); A_joint(2,4); A_joint(4,4);
              A_joint(1,2); A_joint(3,2); A_joint(2,2); A_joint(4,2)];
           
    %solved position of point P  
    p.h1 = Z_leng(1);   
    p.d1 = Z_leng(2);
    
    disp('Initial conditions created, solving EOM...');
    
%% SOLVING EOM AT EACH TIME STEP []
    %creating time array
    tfinal = 5;            % the final time value [s]
    deltat = .0001;        % the time step [s]
    t = 0:deltat:tfinal;   % time array for evaluating solution [s]

    %solve the EoM with ode45
    fdynamicDAE = @(t,X) RSI_dynamics_DAE(t,X,p);   %DAE approach
    fdynamicEL = @(t,X) RSI_dynamics_EL(t,X,p);     %EL approach
    options = odeset('absTol',1e-12,'relTol',1e-12);
    
    if(EL == false && DAE == true)
        [t,X] = ode45(fdynamicDAE,t,X0_DAE,options);
    elseif(EL == true && DAE == false)
        [t,X] = ode45(fdynamicEL,t,X0_EL,options);
    else
        disp('error: must choose a calculation method, DAE or EL')
    end

    disp('EOM solved, post-processing results...');

%% KINEMATIC POST-PROCESSING []
    %get other important kinematic quantities from C.O.M. data
    [Z_E, ~, Z_H, ~, Z_I] = kinematic_processing_finals(X,p);

    %extract C.O.M. position, velocity, angle, & angular velocity @ each 
    %time step
    th_ED = X(:,1);     %link ED state
    thdot_ED = X(:,2); 
    x_G1 = X(:,3);
    xdot_G1 = X(:,4);
    y_G1 = X(:,5);
    ydot_G1 = X(:,6);
    th_EF = X(:,7);     %link EF state
    thdot_EF = X(:,8);
    x_G2 = X(:,9);
    xdot_G2 = X(:,10);
    y_G2 = X(:,11);
    ydot_G2 = X(:,12);
    th_HI = X(:,13);    %link HI state
    thdot_HI = X(:,14);
    x_G3 = X(:,15);
    xdot_G3 = X(:,16);
    y_G3 = X(:,17);
    ydot_G3 = X(:,18);
    r_PI = X(:,19);     %end effector state
    rdot_PI = X(:,20);
    x_P = X(:,21);      %point P state  
    xdot_P = X(:,22);
    y_P = X(:,23);
    ydot_P = X(:,24);
    x_F = X(:,25);      %point F state  
    xdot_F = X(:,26);
    y_F = X(:,27);
    ydot_F = X(:,28);
    
    %point E
    x_E = Z_E(:,1);
    y_E = Z_E(:,2);
    xdot_E = Z_E(:,3);
    ydot_E = Z_E(:,4);
    
    %point H
    x_H = Z_H(:,1);
    y_H = Z_H(:,2);
    xdot_H = Z_H(:,3);
    ydot_H = Z_H(:,4);
    
    %point I
    x_I = Z_I(:,1);
    y_I = Z_I(:,2);
    xdot_I = Z_I(:,3);
    ydot_I = Z_I(:,4);
    
%% SYSTEM ANIMATION []
    %setting up the animation frame
    figure(1)
    hold on
    ylim([-0.4 1])
    xlim([-0.6 1])
    plot(0,0,'r.','LineWidth', 1.5);    %origin (point D)                                     

    %abdominal wall
    abd1 = plot([p.d1-0.1, p.d1-0.04], [-p.h1,-p.h1], 'LineWidth', 1.5);    
    abd1.Color = '#805500';
    abd2 = plot([p.d1+.04, 1], [-p.h1,-p.h1], 'LineWidth', 1.5);           
    abd2.Color = '#805500';
    abd3 = plot([p.d1-0.1, p.d1-0.1], [-p.h1,-1], 'LineWidth', 1.5);    
    abd3.Color = '#805500';
    entry = plot(p.d1,-p.h1,':s','MarkerSize',20);      %point of entry
    entry.Color = '#805500';
    
    %finalizing
    xlabel('X-axis (m)');
    ylabel('Y-Axis (m)');
    title('Robotic Surgery Instrument Dynamics')
    
    %animation processing
    for j = 1:(length(t)-1)/20+1
        i = (j-1)*20+1;

        %first linkage               
        h1 = plot([0,x_E(i)],[0,y_E(i)],'k-', 'LineWidth', 4);              %linkage ED 
        h3 = plot(x_E(1:i),y_E(1:i),'c-');                                  %time trace for point E
        
        %second linkage
        h4 = plot([x_E(i),x_F(i)],[y_E(i),y_F(i)],'k-', 'LineWidth', 4);    %linkage EF
        h6 = plot(x_F(1:i),y_F(1:i),'r-');                                  %time trace for point F
        
        %third linkage
        h7 = plot([x_H(i),x_P(i)],[y_H(i),y_P(i)],'k-', 'LineWidth', 4);    %linkage HP
        
        %end effector
        h10 = plot([x_P(i),x_I(i)],[y_P(i),y_I(i)],'m-', 'LineWidth', 2.5); %linkage PI (end effector)
        h11 = plot(x_I(1:i),y_I(1:i),'g-');                                 %time trace for point I
    
        pause(deltat);     %tuning animation speed
        
        if i==length(t)
            break
        end
        
        %replacing each frame with the next
        delete([h1]);        
        delete([h4]);
        delete([h7]);
        delete(h10);
    end

%% ERROR PLOTTING []
    %linkage length error figures
    error_ED = p.el_ED*ones(length(t),1) - sqrt(x_E.^2 + y_E.^2);
    error_EF = p.el_EF*ones(length(t),1) - sqrt((x_F-x_E).^2 + (y_F - y_E).^2);
    error_HF = p.el_HF*ones(length(t),1) - sqrt((x_F-x_H).^2 + (y_F - y_H).^2);
    error_FP = p.el_FP*ones(length(t),1) - sqrt((x_F-x_P).^2 + (y_F - y_P).^2);
    error_HP = el_HP*ones(length(t),1) - sqrt((x_H-x_P).^2 + (y_H - y_P).^2);
    error_P = sqrt(p.h1^2 + p.d1^2)*ones(length(t),1) - sqrt(x_P.^2 + y_P.^2);
    
    figure(2)
    tiledlayout(1,2)
    nexttile
    plot(t, error_ED)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in length of ED')
    nexttile
    plot(t, error_EF)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in length of EF')
    
    figure(3)
    tiledlayout(1,2)
    nexttile
    plot(t, error_P)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in position of point P')
    nexttile
    plot(t, error_HP)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in length of HP')
    
    figure(4)
    tiledlayout(1,2)
    nexttile
    plot(t, error_HF)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in length of HF')
    nexttile
    plot(t, error_FP)
    xlabel('time(s)')
    ylabel('error (m)')
    title('error in length of FP')

%% ADDITIONAL PLOTTING
    %setting up the animation frame
    figure(1)
    hold on
    ylim([-0.4 1])
    xlim([-0.6 1])
    plot(0,0,'r.','LineWidth', 1.5);    %origin (point D)                                     

    %abdominal wall
    abd1 = plot([p.d1-0.1, p.d1-0.04], [-p.h1,-p.h1], 'LineWidth', 1.5);    
    abd1.Color = '#805500';
    abd2 = plot([p.d1+.04, 1], [-p.h1,-p.h1], 'LineWidth', 1.5);           
    abd2.Color = '#805500';
    abd3 = plot([p.d1-0.1, p.d1-0.1], [-p.h1,-1], 'LineWidth', 1.5);    
    abd3.Color = '#805500';
    entry = plot(p.d1,-p.h1,':s','MarkerSize',20);      %point of entry
    entry.Color = '#805500';
    
    %finalizing
    xlabel('X-axis (m)');
    ylabel('Y-Axis (m)');
    title('Robotic Surgery Instrument Dynamics')

    i = length(t);
    
    %first linkage               
    h1 = plot([0,x_E(i)],[0,y_E(i)],'k-', 'LineWidth', 4);              %linkage ED
%         h2 = plot(x_G1(i),y_G1(i),'wx','MarkerSize',5);                     %C.O.M. 1   
    h3 = plot(x_E(1:i),y_E(1:i),'c-');                                  %time trace for point E

    %second linkage
    h4 = plot([x_E(i),x_F(i)],[y_E(i),y_F(i)],'k-', 'LineWidth', 4);    %linkage EF
%         h5 = plot(x_G2(i),y_G2(i),'wx','MarkerSize',5);                     %C.O.M. 2
    h6 = plot(x_F(1:i),y_F(1:i),'r-');                                  %time trace for point F

    %third linkage
%         h7 = plot([x_H(i),x_F(i)],[y_H(i),y_F(i)],'k-', 'LineWidth', 4);  %linkage HP
    h7 = plot([x_H(i),x_P(i)],[y_H(i),y_P(i)],'k-', 'LineWidth', 4);    %linkage HP
%         h9 = plot(x_G3(i),y_G3(i),'wx','MarkerSize',5);                     %C.O.M. 3

    %end effector
    h10 = plot([x_P(i),x_I(i)],[y_P(i),y_I(i)],'m-', 'LineWidth', 2.5); %linkage PI (end effector)
    h11 = plot(x_I(1:i),y_I(1:i),'g-');                                 %time trace for point I
   
    
