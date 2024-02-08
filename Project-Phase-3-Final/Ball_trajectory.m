%%%%%%%%%%%%%%%%%%%%%%%%%% ME EN 2030 Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mfile to calculate the trajectory of a baseball given environmental %%%
%%% parameters, ball parameters, and initial trajectory values from     %%%
%%% measurements or the solution of the corresponding ball-bat impact   %%%
%%% problem.                                                            %%%
%%%%%%%%%%%%%%%%%%%%%%% Version date: 4/16/2023 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear;clc
%%% Define parameters typical specified units are used for most %%%
% environmental
rho_a = 0.0752; %[lb/ft^3] air density at sea level for Inter Standard Atmosphere
g     = 32.2;  %[ft/sec^2] gravity

% ball and bat
m_ball = 5; %[oz] mass of a baseball (5-5.25 oz for pro baseball)
M_Bat = 30; %[oz] mass of the bat (typically a few oz less than the length)
C = 9; %[in] circumference (9-9.25 in for a pro baseball)
D_Bat = 2.75; %[in] diameter of the baseball bat at impact point

% tee variables
x_o  = 0; %[ft] initial position
y_o  = 3; %[ft] initial height (i.e., the tee height)

%%%%% These parameters are leftover from the original implementation from
%%%%% Phase I, they can be uncommented (and the ball-bat part commented to
%%%%% go back to Phase I initial testing.
%V_bo = [20:0.5:60]; %[mph] initial ball velocity range as [min:step:max]
%theta_o = [10 20 30 40 50]; %[degrees] ball launch angle (relative to horizontal)
%omega = 1800; %[rpm] ball angular velocity (backspin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters from ball-bat impact problem (replaces above)
V_Bat = 30.9; %[mph] bat speed at impact 
theta_Bat = 40; %[degrees] bat angle from horizontal
offset = 1.25; %[in] vertical offset between ball and bat centers of mass at impact
epsilon = 0.69; %[unitless] coefficient of restitution between ball and bat
mu = 0.1; %[unitless] coefficient of friction between ball and bat

%check that offset isn't too large
if(offset >= D_Bat/2+C/(2*pi))
    disp('-------------------------------------------------')
    disp('Oh that is a swing and miss! Try different offset')
    disp('-------------------------------------------------')
    return
end
%initialize variables
V_bo = zeros(length(V_Bat),length(theta_Bat));
theta_o = zeros(length(V_Bat),length(theta_Bat));
omega = zeros(length(V_Bat),length(theta_Bat));
%loop through values
for m=1:length(theta_Bat)
    for n=1:length(V_Bat)
        [V_bo(n,m),theta_o(n,m),~,~,omega(n,m),~] = get_batballimpact(...
            V_Bat(n),theta_Bat(m),offset,epsilon,D_Bat,M_Bat,C,m_ball,mu);
    end
end

%% Integration of the flight of the ball
%initialize variables
max_x = zeros(length(V_Bat),length(theta_Bat));
%loop through values
for m=1:length(theta_Bat)
    for n=1:length(V_Bat)
        [x_b,y_b,Vball] = get_ballflight(x_o,y_o,V_bo(n,m),theta_o(n,m),...
            omega(n,m),C,m_ball,rho_a,g);
        max_x(n,m) = x_b(end);
    end
end
%% Basic Plot of Distance and Speed For Max Speed Tested%%
figure;
%plot the trajectory on the left axis
yyaxis left
plot(x_b,y_b,'LineWidth',2)
ylabel('Height [ft]')
YL = get(gca,'Ylim');YL(1)=0;
set(gca,'Ylim',YL,'FontSize',16)
%plot the baseball's speed (converted to mph) on the right
yyaxis right
plot(x_b,Vball*3600/5280,'LineWidth',2)
ylabel('Speed [mph]')
xlabel('Distance [ft]')
set(gca,'Xlim',[0 max(x_b)],'FontSize',16)
title(['Flight trajectory and speed for V_{Bat}= ',num2str(V_Bat(end)),' mph']);grid on

figure;hold on
LG={};
for m=1:length(theta_Bat)
    plot(max_x(:,m),V_Bat,'-o','LineWidth',2,'MarkerSize',8)
    LG = cat(1,LG,['\theta_{bat} = ',num2str(theta_Bat(m)),'^o']);
end
ylabel('Bat speed [mph]');xlabel('Distance in air [ft]')
title('Distance vs bat speed for different \theta_{bat} values');
legend(LG,'Location','NorthWest')
set(gca,'FontSize',16);grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Function to calculate the ball/bat collision %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V_bo,theta_o,V_Bat_AF,theta_Bat_AF,omega_ball,omega_Bat] = ...
    get_batballimpact(V_Bat,theta_Bat,offset,epsilon,D_Bat,M_Bat,C,m_ball,mu)

%convert all parameters to common units
V_Bat = V_Bat*5280/3600;
theta_Bat = deg2rad(theta_Bat);
offset = offset/12;
R_Bat = D_Bat/2/12;
r_ball = C/(2*pi)/12;

% calculate some initial variables
psi = asin(offset/(R_Bat+r_ball)); %angle to impact plane
V_Batx = V_Bat*cos(theta_Bat);     %x-component of bat velocity
V_Baty = V_Bat*sin(theta_Bat);     %y-component of bat velocity
%solve system of equations

%Coefficient Matrix:
A = [m_ball/M_Bat 0            1         0         0          0; ...
     0            m_ball/M_Bat 0         1         0          0; ...
     -sin(psi)    cos(psi)     0         0         2/5*r_ball 0; ...
     0            0            sin(psi)  -cos(psi) 0          1/2*R_Bat; ...
     cos(psi)     sin(psi)     -cos(psi) -sin(psi) 0          0; ...
     -sin(psi)    cos(psi)     sin(psi)  -cos(psi) -r_ball    -R_Bat];

%Forcing vector
b = [V_Batx;...
     V_Baty;...
     0;...
    -V_Baty*cos(psi)+V_Batx*sin(psi);...
     epsilon*(V_Batx*cos(psi)+V_Baty*sin(psi));...
     0];

%solution vector using Matlab linear solve ('slash' command)
x = A\b;

%check for friction
if(mu*(x(1)*cos(psi)+x(2)*sin(psi)) < x(1)*sin(psi)-x(2)*cos(psi))
    disp('using friction')
    %Modify coefficient matrix:
    A(6,:) = [mu*(cos(psi)-sin(psi)) mu*sin(psi)+cos(psi) 0 0 0 0];

    %re-solve using the new coefficient matrix
    x = A\b;
end

%set outputs with desired units
V_bo = sqrt(x(1)^2+x(2)^2)*3600/5280; %[mph]
theta_o = rad2deg(atan(x(2)/x(1))); %[degrees]
V_Bat_AF = sqrt(x(3)^2+x(4)^2)*3600/5280; %[mph]
theta_Bat_AF = rad2deg(atan(x(4)/x(3))); %[degrees]
omega_ball = x(5)*60/(2*pi); %[rpm]
omega_Bat = x(6)*60/(2*pi); %[rpm]

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Function to integrate the trajectory %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_b,y_b,Vball] = get_ballflight(x_o,y_o,V_bo,theta_o,omega,C,m_b,rho_a,g)

% numerical parameters
dt = 0.01; %[sec] integration timestep (smaller values increase numerical accuracy)

%precomputed variables (can be done before loop)
A_b = 1/(4*pi)*(C/12)^2; %[ft] cross sectional area
C = C/12; %circumference converted to ft
m_b = m_b/16; %ball mass converted to lbs
omega = omega/60; %angular velocity in rev/sec

%initialize variables for loop
y_b = y_o; 
x_b = x_o;
Vball = V_bo*(5280/3600); %converted to ft/sec
u_b = Vball*cosd(theta_o);
v_b = Vball*sind(theta_o);

while y_b(end) > 0 %until the ball hits the ground
    %advance position using last velocity
    x_b = [x_b,x_b(end)+u_b*dt]; %#ok<AGROW> 
    y_b = [y_b,y_b(end)+v_b*dt]; %#ok<AGROW> 

    Vball = [Vball,sqrt(u_b^2+v_b^2)]; %#ok<AGROW> 
    CD = get_CD(Vball(end)); %drag coefficient
    CL = get_CL(Vball(end),C,omega); %lift coefficient
    V_rlx =  (rho_a*A_b*Vball(end))/(2*m_b);

    %advance velocity
    u_b = u_b-V_rlx*(CD*u_b+CL*v_b)*dt;
    v_b = v_b+(V_rlx*(CL*u_b-CD*v_b)-g)*dt;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lift and Drag coefficient functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Drag coefficient function based on Adair 2002 ('Physics of Baseball, 3rd Ed.) 
function CD = get_CD(Vin)
Vb = Vin*3600/5280; %convert to mph for fit
if(Vb < 40)
    CD = 0.5;
else
    CD = -5.981763e-9*Vb^4+2.47429e-6*Vb^3 - ...
        3.383981e-4*Vb^2+1.521796e-2*Vb+0.2917233;
end
end

%Lift coefficient function based on 
function CL = get_CL(Vin,C,omega)
S = C*omega/Vin;
if(S < 0.1)%model of SHS for a baseball
    CL = 1.5*S;
else
    CL = 0.09+0.6*S;
end
end