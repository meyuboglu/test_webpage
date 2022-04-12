%% MPC project Mert Eyuboglu - 1599925
clc, clear all,

%% -------------------------- Model Parameters ------------------------- %%

 R = 1;         %   [Ohm] Electrical resistance
 L = 1/1000;    %   [H] Electrical inductance
 k = 6/100;     %   [NA^-1] Motor constant
 b = 1/1000;    %   [Nsm^-1] Friction coefficent
 m = 7/100;     %   [kg] Pendulum mass 
 l = 1/10;      %   [m] Pendulum length
 J = m*l^2;     %   [kgm^2] Pendulum inertia
 g = 9.81;      %   [ms^-2] Standard gravity
%% ------------------------- 2.a) Linearization --------------------------%%
syms x1 x2 x3 u
x_ss = [1.5707; 0; -1.1445]; u_ss = -1.1445;

% system dynamics as differential equations
f1 = x2;
f2 = (m*g*l/J)*sin(x1)-(b/J)*x2+(k/J)*x3;
f3 = -(k/L)*x2-(R/L)*x3+(1/L)*u;

% Jacobians
J_x = jacobian([f1, f2, f3], [x1, x2, x3]);
J_u = jacobian([f1, f2, f3], u);

% linearized ss model in ct
x1 = x_ss(1); x2 = x_ss(2); x3 = x_ss(3);
u = u_ss;
Ac = double(subs(J_x));
Bc = double(subs(J_u));
Cc = eye(3);
Dc = 0;

% dt ss model
sysc = ss(Ac,Bc,Cc,Dc);
Ts = 0.004;
sysd = c2d(sysc,Ts);
A = sysd.A;
B = sysd.B;
C = sysd.C;
D = sysd.D;

% affine terms
E = [A-eye(3) B; C zeros(3,1)]*[x_ss; u_ss];
g1 = -E(1:3);
g2 = x_ss - E(4:end);

%% -------------------------- 2.b-c) Linear MPC -------------------------%%
% cost function
Q = [5000 0 0; 0 1 0; 0 0 1]; %change for tuning
R = 0; % change for tuning
[K,P_lqr,~] = dlqr(A,B,Q,R);
K_lqr = -K;
N = 4; % prediction horizon

% prediction matrices
[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P_lqr,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega*Phi;

% constrained MPC
xmax = [2*pi; 18; 6];
xmin = [-2*pi; -18; -6];
umax = 8;
umin = -8;

% closed loop lqr
A_cl = A + B*K_lqr;

% Compute the invariant set with constraint admissible set
m = size(B,2);
n = size(B,1);
model = LTISystem('A', A_cl,'Ts',Ts);
% state constraints (A_x * x <= b_x)
X_set = Polyhedron([eye(n);-eye(n)],[xmax; -xmin])-x_ss;
% input constraints (A_u * u <= b_u)
U_set = Polyhedron([eye(m); -eye(m)],[umax; -umin])-u_ss;
% constraints admissible set
CA_set = Polyhedron([eye(m);-eye(m)]*K_lqr,[umax; -umin]);
% input constraint admissible set
IA_set = X_set&CA_set;
% invariant set with CA set
INVset_mpc = model.invariantSet('X',IA_set);

% feasible set of states
MN_mpc = INVset_mpc.A;
bN_mpc = INVset_mpc.b;

[W, LL, c] = getWLc_TS( A, B, xmax-x_ss, xmin-x_ss, umax-u_ss, umin-u_ss, Gamma, Phi, MN_mpc, bN_mpc);

% projection of feasible set of states on x1 and x2, also x0 marked
x0 = [5.9832; 17.7; -1.549];
P = Polyhedron('A',[-W LL],'B',c);
Feasiblesetx = projection(P,1:3);
figure()
plot(x0(1),x0(2),'k*')
hold on
plot(projection(Feasiblesetx+x_ss,1:2))
title('Projection of the feasible set of states on x1 and x2')
xlabel('x1') 
ylabel('x2')

% try a feasible point to see if the controller is working 
x0 = [1.8708; 0; -1.549];
P = Polyhedron('A',[-W LL],'B',c);
Feasiblesetx = projection(P,1:3);
figure()
plot(projection(Feasiblesetx+x_ss,1:2))
hold on
plot(x0(1),x0(2),'k*')

% Trajectory constrained MPC
k_sim = 500;
x = zeros(size(B,1),k_sim+1);
x(:,1) = x0;
options_qp =  optimoptions('quadprog','Display','off');
for i = 1:k_sim
    x_k = x(:,i);
    z_k = x_k - x_ss;
    [u_qp,fval,exitflag] = quadprog(G,F*z_k,LL,c+W*z_k,[],[],[],[],[],options_qp);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exit == -2
            fprint('Optimization problem is infeasible. \n')
        end
    end
    u_mpc(i) = u_qp(1)+u_ss;
    x(:,i+1) = A*x(:,i) + B*u_mpc(i) + g1;
end
figure;
subplot(2,2,1);
stairs(0:k_sim,x(1,:))
hold on
stairs(0:k_sim,x(2,:))
hold on
stairs(0:k_sim,x(3,:))
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
title('Constrained Linear MPC');
legend({'x1(angle)','x2(velocity)','x3(current)'},'Location','northeast')
subplot(2,2,2);
stairs(0:k_sim-1,u_mpc)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('Constrained Linear MPC');

%% simulate linear MPC with the nonlinear Simulink model
clear u; clear x;

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k = 6/100;     %   [NA^-1] Motor constant
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity

T_sample = 0.004;
T_sim = 1;
k_sim = 500;
x = zeros(size(B,1),k_sim+1);
x(:,1) = x0;
options_qp =  optimoptions('quadprog','Display','off');
for i = 1:k_sim
    tic
    x_k = x(:,i);
    z_k = x_k - x_ss;
    [u_qp,fval,exitflag] = quadprog(G,F*z_k,LL,c+W*z_k,[],[],[],[],[],options_qp);
    if exitflag ~= 1
        warning('exitflag quadprog =%d\n', exitflag)
        if exit == -2
            fprint('Optimization problem is infeasible. \n')
        end
    end
    u_mpc(i) = u_qp(1)+u_ss;
    u_input = u_mpc(i);
    x1_0 = x(1,i);
    x2_0 = x(2,i);
    x3_0 = x(3,i);
    time_lin(i) = toc;
    sim('Pendulum_Nonlinear_System')
    x(:,i+1) = [x1_state(end); x2_state(end); x3_state(end)];
end

figure;
subplot(2,2,1);
stairs(0:k_sim,x(1,:))
hold on
stairs(0:k_sim,x(2,:))
hold on
stairs(0:k_sim,x(3,:))
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
title('Constrained Linear MPC with nonlinear model');
legend({'x1(angle)','x2(velocity)','x3(current)'},'Location','northeast')
subplot(2,2,2);
stairs(0:k_sim-1,u_mpc)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('Constrained Linear MPC with nonlinear model');

% average computation toime
ave_time_linear = sum(time_lin)/k_sim; 

%% bonus question
for i = 1:5
    N = i;
    [Phi, Gamma] = ABN2PhiGamma(A,B,N);
    [ W, L, c] = getWLc_TS( A, B, xmax-x_ss, xmin-x_ss, umax-u_ss, umin-u_ss, Gamma, Phi, MN_mpc, bN_mpc);
    P = Polyhedron('A',[-W L],'B',c);
    Feasiblesetx(i) = projection(P,1:3);
    plot(projection(Feasiblesetx(i)+x_ss,1:2),'color','yellow','Alpha',1-0.9)
    hold on
end
plot(5.9832,17.7,'k*')

%% lpv nonlin mpc
clear x_ss; clear u_ss; clear g1; clear g2; clear x; clear u_mpc;

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k = 6/100;     %   [NA^-1] Motor constant
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity

N=4;
Ts = 0.004;
A = [1 Ts 0; 0 1-Ts*(b/J) Ts*(k/J); 0 -Ts*(k/L) 1-Ts*(R/L)];
B = [0; 0; Ts/L];
C = [1 0 0];
% affine term: g = [0; Ts*(m*g*l/J)*sin(x1); 0]
x0 = [4.7124; 17.7; -1.549];

[Phi, Gamma] = ABN2PhiGamma(A,B,N);
[Psi, Omega] = QRPN2PsiOmega(Q,R,P_lqr,N);
G = 2*(Psi+Gamma'*Omega*Gamma);
F = 2*Gamma'*Omega;

xmax = [2*pi; 18; 6];
xmin = [-2*pi; -18; -6];
umax = 8;
umin = -8;
[D, M, E, c] = getDMEc( A, B, xmax, xmin, umax, umin, Gamma, Phi);
LL = M*Gamma+E;
W = -D-M*Phi;

% Trajectory constrained MPC
k_sim =500;
x = zeros(3,k_sim+1);
x(:,1) = x0;
r = 1.5707;
X = [x0 x0 x0 x0];%initial estimated trajectory
options_qp =  optimoptions('quadprog','Display','off');
for i = 1:k_sim
    % iterate at each time step
    for iter = 1:Inf
        % compute affine term x_ss, u_ss
        g1(:,1) = [0; Ts*(m*g*l/J)*sin(x(1,i)); 0];
        x_ss_u_ss(:,1) = (inv([A-eye(3),B;C,0]))*[-g1(:,1);r];
        for j = 2:N
            g1(:,j) = [0; Ts*(m*g*l/J)*sin(X(1,j-1)); 0];
            x_ss_u_ss(:,j) = (inv([A-eye(3),B;C,0]))*[-g1(:,j);r];
        end
        u_ss = x_ss_u_ss(4,:)';
        x_ss = [x_ss_u_ss(1:3,1);x_ss_u_ss(1:3,2);x_ss_u_ss(1:3,3);x_ss_u_ss(1:3,4)];
        
        % compute prediciton matrix corresponding to the affine term
        T = g1(:,1);
        for ii = 2:1:N
            T = [T; A*T(end-2:end,:) + g1(:,ii)];
        end
        T = double(T);
        
        % solve qp
        [u_qp,fval,exitflag] = quadprog(G,F*Phi*x(:,i)+F*T-F*x_ss-2*Psi*u_ss,LL,c+W*x(:,i)-M*T,[],[],[],[],[],options_qp);
        if exitflag ~= 1
            warning('exitflag quadprog =%d\n', exitflag)
            if exit == -2
                fprint('Optimization problem is infeasible. \n')
            end
        end
        
        % update future states by applying predicted input
        X(:,1) = A*x(:,i) + B*u_qp(1) + g1(:,1);
        for n = 2:N
            X(:,n) = A*X(:,n-1) + B*u_qp(n) + [0; Ts*(m*g*l/J)*sin(X(1,n-1)); 0];
        end
        u_qp_list(:,iter) = u_qp;
        % check if parameters estimated accurately
        if iter>1 && norm(u_qp_list(:,iter)-u_qp_list(:,iter-1),Inf)<1e-5
            break
        end
    end
    u_mpc(i) = u_qp(1);
    x(:,i+1) = A*x(:,i) + B*u_mpc(i)+ g1(:,1);
    X = double(X(:,2:end));
end

figure;
subplot(2,2,1);
stairs(0:k_sim,x(1,:))
hold on
stairs(0:k_sim,x(2,:))
hold on
stairs(0:k_sim,x(3,:))
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
title('constrained MPC');
subplot(2,2,2);
stairs(0:k_sim-1,u_mpc)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('constrained MPC');

%% simulate LPV MPC with nonlinear Simulink model
clear x_ss; clear u_ss; clear g1; clear g2; clear x; clear u_mpc;

R = 1;         %   [Ohm] Electrical resistance
L = 1/1000;    %   [H] Electrical inductance
k = 6/100;     %   [NA^-1] Motor constant
b = 1/1000;    %   [Nsm^-1] Friction coefficent
m = 7/100;     %   [kg] Pendulum mass
l = 1/10;      %   [m] Pendulum length
J = m*l^2;     %   [kgm^2] Pendulum inertia
g = 9.81;      %   [ms^-2] Standard gravity

% Trajectory constrained MPC
T_sample = 0.004;
T_sim = 1;
k_sim =500;
x = zeros(3,k_sim+1);
x(:,1) = x0;
r = 1.5707;
X = [x0 x0 x0 x0];
options_qp =  optimoptions('quadprog','Display','off');
for i = 1:k_sim
    tic
    for iter = 1:Inf
        g1(:,1) = [0; Ts*(m*g*l/J)*sin(x(1,i)); 0];
        x_ss_u_ss(:,1) = (inv([A-eye(3),B;C,0]))*[-g1(:,1);r];
        for j = 2:N
            g1(:,j) = [0; Ts*(m*g*l/J)*sin(X(1,j-1)); 0];
            x_ss_u_ss(:,j) = (inv([A-eye(3),B;C,0]))*[-g1(:,j);r];
        end
        u_ss = x_ss_u_ss(4,:)';
        x_ss = [x_ss_u_ss(1:3,1);x_ss_u_ss(1:3,2);x_ss_u_ss(1:3,3);x_ss_u_ss(1:3,4)];

        T = g1(:,1);
        for ii = 2:1:N
            T = [T; A*T(end-2:end,:) + g1(:,ii)];
        end
        T = double(T);

        [u_qp,fval,exitflag] = quadprog(G,F*Phi*x(:,i)+F*T-F*x_ss-2*Psi*u_ss,LL,c+W*x(:,i)-M*T,[],[],[],[],[],options_qp);
        if exitflag ~= 1
            warning('exitflag quadprog =%d\n', exitflag)
            if exit == -2
                fprint('Optimization problem is infeasible. \n')
            end
        end
        
        X(:,1) = A*x(:,i) + B*u_qp(1) + g1(:,1);
        for n = 2:N
            X(:,n) = A*X(:,n-1) + B*u_qp(n) + [0; Ts*(m*g*l/J)*sin(X(1,n-1)); 0];
        end
        u_qp_list(:,iter) = u_qp;
        if iter>1 && norm(u_qp_list(:,iter)-u_qp_list(:,iter-1),Inf)<1e-5
            break
        end
    end
    u_mpc(i) = u_qp(1);
    u_input = u_mpc(i);
    x1_0 = x(1,i);
    x2_0 = x(2,i);
    x3_0 = x(3,i);
    time_lpv(i) = toc;
    sim('Pendulum_Nonlinear_System')
    x(:,i+1) = [x1_state(end); x2_state(end); x3_state(end)];
    X = double(X(:,2:end));
end

figure;
subplot(2,2,1);
stairs(0:k_sim,x(1,:))
hold on
stairs(0:k_sim,x(2,:))
hold on
stairs(0:k_sim,x(3,:))
xlabel('$k$','Interpreter','latex');
ylabel('$x$','Interpreter','latex');
title('LPV MPC with nonlinear Simulink model');
legend({'x1(angle)','x2(velocity)','x3(current)'},'Location','northeast')
subplot(2,2,2);
stairs(0:k_sim-1,u_mpc)
xlabel('$k$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('LPV MPC with nonlinear Simulink model');

% average computation time
ave_time_lpv = sum(time_lpv)/k_sim; 








