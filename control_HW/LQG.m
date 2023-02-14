syms  g M m_1 m_2 l_2 l_1 
%Constants

A=[0 1 0 0 0 0;
   0 0 (-g*m_1)/M 0 (-g*m_2)/M 0;
   0 0 0 1 0 0;
   0 0 (-g*(M+m_1))/(l_1*M) 0 (-g*m_2)/(M*l_1) 0;
   0 0 0 0 0 1;
   0 0 (-g*m_1)/(M*l_2) 0 (-g*(M+m_2))/(M*l_2) 0];
B=[0;1/M;0;1/(M*l_1);0;1/(M*l_2)];
D=0;
C_1=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

Q = eye(6,6);
Q(1,1)=10;
Q(2,2)=1;
Q(3,3)=100;
Q(4,4)=1000;
Q(5,5)=1;
Q(6,6)=800;
R= 0.0001;
X_0 = [2; 0; 0.1; 0; 0.1; 0];

%Substituting
A_1=double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
B_1 = double(subs(B,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));

[K1,S,~] = lqr(A_1, B_1, Q, R);

sys = ss(A_1-B_1*K1,B_1,C_1,D);
X_final = [0;0;0;0;0;0];

%Kalman Estimator
Bd = 0.01*eye(6); %disturbance
Vn = 0.001;     %Gaussian White Noise
U = @(state) -K1*(X_0 - X_final);
[L,P,E] = lqe(A_1,Bd,C_1,Bd,Vn*eye(3)); 
Ac1 = A_1-(L*C_1);
e_sys1 = ss(Ac1,[B_1 L],C_1,0);
ts = 0:0.01:100;
[t,state1] = ode45(@(t,state)nonLinear1_LQG(t,state,U(state),L),ts,X_0);
figure();
hold on
plot(t,state1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear LQG for output vector: x(t)')
legend('x')
hold off




%% LQG Non Linear
function ds1 = nonLinear1_LQG(t,x,f,l1)
    m1 = 100; m2 = 100; M=1000; L1 = 20; L2 = 10; g = 9.81;
    X_ = x(1);
    X_dot = x(2);
    theta_1 = x(3);
    theta_dot_1 = x(4);
    theta_2 = x(5);
    theta_dot2 = x(6);
    ds1 = zeros(6,1);
    y1 = [X_;0;0];
    c_1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
    iter = l1*(y1-c_1*x);
    ds1(1) = X_dot+iter(1);
    ds1(2) = (f - m1*L1*sin(theta_1)*theta_dot_1^2 - (m2*L2*sin(theta_2)*theta_dot2^2 - m1*g*sin(theta_1)*cos(theta_1) - m2*g*sin(theta_2)*cos(theta_2))/(M + m1 + m2 - m1*cos(theta_1)^2 - m2*cos(theta_2)^2))+iter(2);
    ds1(3) = theta_dot_1+iter(3);
    ds1(4) = cos(theta_1)*ds1(2)/L1 - (g*sin(theta_1)/L1)+iter(4);
    ds1(5) = theta_dot2+iter(5);
    ds1(6) = cos(theta_2)*ds1(2)/L2 - (g*sin(theta_2)/L2)+iter(6);
end