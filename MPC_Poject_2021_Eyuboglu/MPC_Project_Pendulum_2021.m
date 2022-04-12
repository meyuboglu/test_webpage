% Electrically driven inverted pendulum
% Based on Section 8-5-2 pag 171 from Hanema J., Anticipative MPC for LPV,
% 2018.
% C:\Users\20205236\OneDrive - TU Eindhoven\Documents\MATLAB\MPC_Coure
% ---------------------------------------------------------------------
% Inthis numerical example is to control the angle q of the electrically
% driven inverted pendulum. Define the state vector
% x = [q dq i] = [ x1 x2 x3] where :
%               q ->  angle of the pendulum   [rad]
%               dq -> angular velocity        [rad/s]  
%               di -> motor current           [A]  
% --------------------------------------------------------------------
clc, clear all,
%% ------------------------ Simulation Parameters ---------------------- %%
%           Time
T_step =  1000;     % [s] 
pu = 0.5;           % [fr] step amplitude
T_sample = 0.004;   % [s]
T_sim = 10000/4;      % []

%% -------------------------- Model Parameters ------------------------- %%

 R = 1;         %   [Ohm] Electrical resistance
 L = 1/1000;    %   [H] Electrical inductance
 k = 6/100;     %   [NA^-1] Motor constant
 b = 1/1000;    %   [Nsm^-1] Friction coefficent
 m = 7/100;     %   [kg] Pendulum mass 
 l = 1/10;      %   [m] Pendulum length
 J = m*l^2;     %   [kgm^2] Pendulum inertia
 g = 9.81;      %   [ms^-2] Standard gravity
 
%% ------------------------ Initial States ------------------------- %%

% For assignment 1. change the initial states in this section.


x1_0 = 1e-5;           % [rad]         q       Angle of the pendulum 

x2_0 = 0;               % [rad/s]       dq      Angular Velocity 

x3_0 = 0;               % [A]           i       Motor Current

u_input = 0;            % [V]           u       DC Voltage

%% ---------------------------- Simulation ----------------------------- %%


                    sim('Pendulum_Nonlinear_System')

%% ------------------------------- Plots ------------------------------- %%
t = linspace(0,10,T_sim+1);
figure()
plot(t,x1_state,'-.',t,x2_state,'-',t,x3_state,'--','LineWidth',1.25)
title('State Trajectories for x(0)=[1e-5  0  0]^T')
xlabel('time (seconds)') 
ylabel('amplitude')
legend({'x1(angle)','x2(velocity)','x3(current)'},'Location','northeast')
