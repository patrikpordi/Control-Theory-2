%% Linearization
%State equation for the Jacobian
%Constants
syms F m_1 m_2 M l_1 l_2 g
%X should contain x,theta_1,theta_2,x_dot, theta_1_dot, theta_2_dot
syms x x_dot theta_1 theta_1_dot theta_2 theta_2_dot
%The six functions are the following based on the state variables and the
%equation of motion
f_1 = x_dot;
f_2 = (F-m_1*l_1*sin(theta_1)*theta_1_dot^2- m_2*l_2*sin(theta_2)*theta_2^2 ...
-m_1*g*sin(theta_1)*cos(theta_1)-m_2*g*sin(theta_2)*cos(theta_2))...
/(M+m_1+m_2-m_1*cos(theta_1)^2-m_2*cos(theta_2)^2);
f_3 =theta_1_dot;
f_4 = -g*sin(theta_1)/l_1+cos(theta_1)*f_2/l_1;
f_5 = theta_2_dot;
f_6 = -g*sin(theta_2)/l_2+cos(theta_2)*f_2/l_2;

%The Jacobians and their linearized form using the equilibrium point is the
%following
%A
J_A = jacobian([f_1,f_2,f_3,f_4,f_5,f_6],[x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
A = simplify((subs(J_A, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
    [0 0 0 0 0 0])));
%B
J_B = jacobian([f_1,f_2,f_3,f_4,f_5,f_6],F);
B = simplify((subs(J_B, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
    [0 0 0 0 0 0])));
% The output vector [x, theta_1, theta_2] was selected
%C
J_C = jacobian ([x, theta_1, theta_2], [x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
C = simplify((subs(J_C, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
    [0 0 0 0 0 0])));
%D
J_D = jacobian([x, theta_1, theta_2],F);
D = simplify((subs(J_D, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
    [0 0 0 0 0 0])));

%% Controllability
%Building the controllability matrix from B-A^5*B
Cont = [B, A*B, A^2*B,A^3*B,A^4*B,A^5*B];
%Matrix in Latex format
latex_table = latex(sym(simplify(Cont)));

rank_condition_controllability=simplify(det(Cont));
disp(rank_condition_controllability);
%Condition in Latex format
latex_condition = latex(sym(simplify(det(Cont))));

%% LQR

%Q & R has to be iterated to find the optimal values
Q = eye(6,6);

A_1 = double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
B_1 = double(subs(B,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
C_1 = double(subs(C,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
D_1 = double(subs(D,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));

R = 0.00001;%add
disp(rank(ctrb(A_1,B_1)));%add
Q(1,1)=10;
Q(2,2)=1;
Q(3,3)=100;
Q(4,4)=1000;
Q(5,5)=1;
Q(6,6)=800;
latex_ku = latex(sym(Q));%add
K_1 = lqr(A_1, B_1, Q, R);%add
I_1=eig(A_1-B_1*K_1);%add

X_0=[2;0;0.1;0;-0.1;0];
sys=ss((A_1-B_1*K_1),B_1,C_1,D_1);
initial(sys,X_0)
figure(1)
%step(sys)
xlabel("t [s]")
ylabel("\theta_2                        \theta_1                            x")

%[T,D]=eig(A_1-K_1*B_1);

%disp(diag(real(D)));
%disp(T(:,1));

%% Non_Linear LQR
tspan = 0:0.1:100;
[t,q] = ode45(@(t,q)LQRonNonLinear(t,q,-K_1*q),tspan,X_0);
figure(2);
hold on
plot(t,q(:,1))
plot(t,q(:,3))
plot(t,q(:,5))
ylabel('"\theta_2                        \theta_1                            x"')
xlabel('t [s]')
title('Non-Linear system using LQR controller')
legend('x','\theta_1','\theta_2')
%% Observability

%Output x, C_1=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]
%J_C_1 = jacobian ([x, 0, 0], [x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
%C_1 = simplify((subs(J_C_1, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
%    [0 0 0 0 0 0])));
C_1=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
O_1=[C_1; C_1*A; C_1*A^2; C_1*A^3; C_1*A^4; C_1*A^5];
O_1_rank=rank(O_1);
%Output theta_1, theta_2, C2=[0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]
J_C_2 = jacobian ([0, theta_1, theta_2], [x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
C_2 = simplify((subs(J_C_2, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
    [0 0 0 0 0 0])));

O_2=[C_2; C_2*A; C_2*A^2; C_2*A^3; C_2*A^4; C_2*A^5];
O_2_rank=rank(O_2);
%rank!=6 so not observable

%Output x, theta_2, C3=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0]
%J_C_3 = jacobian ([x, 0, theta_2], [x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
%C_3 = simplify((subs(J_C_3, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
%    [0 0 0 0 0 0])));
C_3=[1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
O_3=[C_3; C_3*A; C_3*A^2; C_3*A^3; C_3*A^4; C_3*A^5];
O_3_rank=rank(O_3);

%Output x, theta_1,theta_2, C4=[1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]
%J_C_4 = jacobian ([x, theta_1, theta_2], [x, x_dot, theta_1, theta_1_dot, theta_2, theta_2_dot]);
%C_4 = simplify((subs(J_C_4, [x x_dot theta_1 theta_1_dot theta_2 theta_2_dot],...
%    [0 0 0 0 0 0])));

C_4=[1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
O_4=[C_4; C_4*A; C_4*A^2; C_4*A^3; C_4*A^4; C_4*A^5];
O_4_rank=rank(O_4);

% not considering system 2 because its not observable
sys_1 = ss(A_1, B_1, C_1, D_1);
sys_3 = ss(A_1, B_1, C_3, D_1);
sys_4 = ss(A_1, B_1, C_4, D_1);

%% Luenberger Observer
% Kalman Filter
Bd = 0.1 * eye(6,6);
Vn = 0.01 * eye(3,3);

[L_1,~, ~] = lqe(A_1, Bd, C_1, Bd,Vn);
[L_3,~, ~] = lqe(A_1, Bd, C_3, Bd,Vn);
[L_4,~, ~] = lqe(A_1, Bd, C_4, Bd,Vn);

A_C1 = A_1 - L_1*C_1;
A_C3 = A_1 - L_3*C_3;
A_C4 = A_1 - L_4*C_4;

%Taking D=0
e_sys_1 = ss(A_C1, [B_1 L_1], C_1, 0);
e_sys_3 = ss(A_C3, [B_1 L_3], C_3, 0);
e_sys_4 = ss(A_C4, [B_1 L_4], C_4, 0);

%%step input response
t_span = 0:0.1:100;
unit_Step = 0*t_span;
unit_Step(200:length(t_span)) = 1;

[Y_1,~] = lsim(sys_1,unit_Step, t_span);
[X_1,~] = lsim(e_sys_1,[unit_Step;Y_1'],t_span);

[Y_3,~] = lsim(sys_3,unit_Step, t_span);
[X_3,~] = lsim(e_sys_3,[unit_Step;Y_3'],t_span);

[Y_4,~] = lsim(sys_4,unit_Step, t_span);
[X_4,t] = lsim(e_sys_4,[unit_Step;Y_4'],t_span);

%Setting up plot
figure();
hold on
plot(t,Y_1(:,1),'y','Linewidth',2)
plot(t,X_1(:,1),'k','Linewidth',1)
ylabel('State Variables->')
xlabel('time(sec)->')
legend('x(t)','Estimated x(t)')
title('Output vector response at step input (x(t)):')
hold off

figure();
hold on
plot(t,Y_3(:,1),'y','Linewidth',2)
plot(t,Y_3(:,3),'b','Linewidth',2)
plot(t,X_3(:,1),'k','Linewidth',1)
plot(t,X_3(:,3),'m--','Linewidth',1)
ylabel('State Variables->')
xlabel('time(sec)->')
legend('x(t)','theta_2(t)','Estimated x(t)','Estimated theta_2(t)')
title('Output vector response at step input (x(t),theta_2(t)):')
hold off

figure();
hold on
plot(t,Y_4(:,1),'y','Linewidth',2)
plot(t,Y_4(:,2),'b','Linewidth',2)
plot(t,Y_4(:,3),'g','Linewidth',2)
plot(t,Y_4(:,1),'k','Linewidth',1)
plot(t,X_4(:,2),'r','Linewidth',1)
plot(t,X_4(:,3),'m--','Linewidth',1)
ylabel('State Variables->')
xlabel('time(sec)->')
legend('x(t)','theta_1(t)','theta_2(t)','Estimated x(t)','Estimated theta_1(t)','Estimated theta_2(t)')
title('Output vector response at step input (x(t),theta_1(t),theta_2(t)):')
hold off

%Linear system 
[t,x1] = ode45(@(t,x)LO_linear1(t, x ,L_1, A_1, B_1, C_1),t_span,X_0);
figure();
hold on
plot(t,x1(:,1))
ylabel('State variables->')
xlabel('time (sec)->')
title('Linear System Luenberger Observer when output vector is x(t):')
legend('x')
hold off

[t,x3] = ode45(@(t,x)LO_linear3(t, x, L_3, A_1, B_1, C_3),t_span,X_0);
figure();
hold on
plot(t,x3(:,1))
plot(t,x3(:,5))
ylabel('State variables->')
xlabel('time (s)->')
title('Linear System Luenberger Observer when output vector  is (x(t),theta_2(t)):')
legend('x','theta_2')
hold off

[t,x4] = ode45(@(t,x)LO_linear4(t, x, L_4, A_1, B_1, C_4),t_span,X_0);
figure();
hold on
plot(t,x4(:,1))
plot(t,x4(:,3))
plot(t,x4(:,5))
ylabel('State variables->')
xlabel('time (s)->')
title('Linear System Luenberger Observer when output vector is (x(t),theta_1(t),theta_2(t)):')
legend('x','theta_1','theta_2')
hold off

%Non-linear system 
[t,q1] = ode45(@(t,q)LO_NonLinear1(t,q,5,L_1),t_span,X_0);
figure();
hold on
plot(t,q1(:,1))
ylabel('State variables->')
xlabel('Time (sec)->')
title('Non-linear System Luenberger Observer when output vector is x(t):')
legend('x')
hold off

[t,q3] = ode45(@(t,q)LO_NonLinear3(t,q,5,L_3),t_span,X_0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('State variables->')
xlabel('Time (s)->')
title('Non-linear System Luenberger Observer when output vector is (x(t),theta_2(t)):')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)LO_NonLinear4(t,q,5,L_4),t_span,X_0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('State variables->')
xlabel('Time (s)->')
title('Non-linear System Luenberger Observer when output vector is (x(t),theta_1(t),theta_2(t)):')
legend('x','theta_1','theta_2')
hold off


%% Function definition
%LQR non linear
function dQ = LQRonNonLinear(~,y,F)
    m1 = 100; m2 = 100; M=1000; L1 = 20; L2 = 10; g = 9.81;
    dx = y(2);
    t1 = y(3);
    td1 = y(4);
    t2 = y(5);
    td2 = y(6);

    %Nonlinear calculations
    dQ = zeros(6,1);
    dQ(1) = dx;
    dQ(2) = (F - m1*L1*sin(t1)*td1^2 - m2*L2*sin(t2)*td2^2 - m1*g*sin(t1)*cos(t1)...
        - m2*g*sin(t2)*cos(t2))/(M + m1 + m2 - m1*cos(t1)^2 - m2*cos(t2)^2);
    dQ(3) = td1;
    dQ(4) = cos(t1)*dQ(2)/L1 - g*sin(t1)/L1;
    dQ(5) = td2;
    dQ(6) = cos(t2)*dQ(2)/L2 - g*sin(t2)/L2;
end

%L observer linear C_1
function dX = LO_linear1(~,x,L,A,B,C)
 y = [x(1); 0; 0];
 R = 0.00001;
 Q(1,1)=10;
 Q(2,2)=1;
 Q(3,3)=100;
 Q(4,4)=1000;
 Q(5,5)=1;
 Q(6,6)=800;
 K = lqr(A, B, Q, R);
 dX = (A-B*K)*x + L*(y - C*x);
end

%L observer linear C_3
function dX = LO_linear3(~,x,L,A,B,C)
 y = [x(1); 0; x(5)];
 R = 0.00001;%add
 Q(1,1)=10;
 Q(2,2)=1;
 Q(3,3)=100;
 Q(4,4)=1000;
 Q(5,5)=1;
 Q(6,6)=800;
 K = lqr(A, B, Q, R);
 dX = (A-B*K)*x + L*(y - C*x);
end

%L observer linear C_4
function dX = LO_linear4(~,x,L,A,B,C)
 y = [x(1); x(3); x(5)];
 R = 0.00001;
 Q(1,1)=10;
 Q(2,2)=1;
 Q(3,3)=100;
 Q(4,4)=1000;
 Q(5,5)=1;
 Q(6,6)=800;
 K = lqr(A, B, Q, R);
 dX = (A-B*K)*x + L*(y - C*x);
end

%L observer non linear C_1
function dQ = LO_NonLinear1(~,y,F,L_2)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    x_d = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y1 = [x; 0; 0];
    c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
    sum = L_2*(y1-c1*y);
    dQ(1) = x_d + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g...
        - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))...
        /(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

%L observer non linear C_3
function dQ = LO_NonLinear3(~,y,F,L_3)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    x_d = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y3 = [x; 0; t2];
    c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
    sum = L_3*(y3-c3*y);
    dQ(1) = x_d + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g...
        - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))...
        /(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

%L observer non linear C_4
function dQ = LO_NonLinear4(~,y,F,L_4)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    x_d = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y4 = [x; t1; t2];
    c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
    sum = L_4*(y4-c4*y);
    dQ(1) = x_d + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g...
        - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))...
        /(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

