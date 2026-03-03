function var_dot = QuadrotorEOM(t, var, g, m, I, d, km, nu, mu, motor_forces)
% for use by ode45 to simulate the full nonlinear equations of motion where: t is time; var is the 12 x
% 1 aircraft state vector; g is the acceleration due to gravity; m is mass; I is the inertia matrix; d, km,
% nu, and mu are the remaining quadrotor parameters; motor_forces = [f1; f2; f3; f4] is
% the 4 x 1 vector of motor forces, and var_dot is the 12 x 1 derivative of the state vector. Include
% attitude dynamics and kinematics using the Euler angle attitude representation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  t: Time 
%  var: 12x1 state vector   [ x, y, z, phi, theta, psi, u, v, w, p, q, r] 
%  g: accel due to grav 
%  m: Mass
%  I: Inertia Matrix
%  d: 
%  km:
%  nu:
%  mu:
%  Motor Forces: [f1; f2; f3; f4] is the 4 x 1 vector of motor forces
%  var_dot: 12x1 derivative of the state vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
State.x=var(1);
State.y=var(2);
State.z=var(3);
State.phi=var(4);
State.theta=var(5);
State.psi=var(6);
State.u=var(7);
State.v=var(8);
State.w=var(9);
State.p=var(10);
State.q=var(11);
State.r=var(12);
%% Trig structure 
% Psi
Trig.cpsi=cos(State.psi);
Trig.spsi=sin(State.psi);
Trig.tpsi=tan(State.psi);
% Theta
Trig.ctheta=cos(State.theta);
Trig.stheta=sin(State.theta);
Trig.ttheta=tan(State.theta);
% Phi
Trig.cphi=cos(State.phi);
Trig.sphi=sin(State.phi);
Trig.tphi=tan(State.phi);

%% Inetria Struct
In.x=I(1,1);
In.y=I(2,2);
In.z=I(3,3);


%% Control Struct
control=[-1 -1 -1 -1; -d./sqrt(2) -d./sqrt(2) d./sqrt(2) d./sqrt(2); d./sqrt(2) -d./sqrt(2) -d./sqrt(2) d./sqrt(2); km -km km -km].*motor_forces;
Control.Z=sum(control(1,:));
Control.L=sum(control(2,:));
Control.M=sum(control(3,:));
Control.N=sum(control(4,:));

%% Aero Struct
Force= -nu.*norm([State.u State.v State.w]).*[State.u; State.v; State.w];
Moments=-mu.*norm([State.p State.q State.r]).*[State.p; State.q; State.r];
Aero.X=Force(1);
Aero.Y=Force(2);
Aero.Z=Force(3);
Aero.L=Moments(1);
Aero.M=Moments(2);
Aero.N=Moments(3);
%% X_dot, Y_dot, Z_dot
Pos_dot=[Trig.ctheta.*Trig.cpsi, (Trig.sphi.*Trig.stheta.*Trig.cpsi)-(Trig.cphi.*Trig.spsi), (Trig.cphi.*Trig.stheta.*Trig.cpsi)-(Trig.sphi.*Trig.spsi); ...
                     Trig.ctheta.*Trig.spsi, (Trig.sphi.*Trig.stheta.*Trig.spsi)+(Trig.cphi.*Trig.cpsi), (Trig.cphi.*Trig.stheta.*Trig.spsi)-(Trig.sphi.*Trig.cpsi); ...
                     -Trig.stheta, Trig.ctheta.*Trig.sphi, Trig.ctheta.*Trig.cphi] * [State.u;State.v;State.w];

X_dot = Pos_dot(1);
Y_dot = Pos_dot(2);
Z_dot = Pos_dot(3);


%% Phi_dot, Theta_dot, Psi_dot
Angle_dot=[ 1, Trig.sphi.*Trig.ttheta, Trig.cphi.*Trig.ttheta;
            0, Trig.cphi, -Trig.sphi
            0, Trig.sphi./Trig.ctheta, Trig.cphi./Trig.ctheta]*[State.p;State.q;State.r;]; 


Phi_dot = Angle_dot(1);
Theta_dot = Angle_dot(2);
Psi_dot = Angle_dot(3);


%% U_dot, V_dot, W_dot

v_dot= cross([State.u,State.v,State.w],[State.p,State.q,State.r])' + g.*[-Trig.stheta;Trig.ctheta.*Trig.sphi ;Trig.ctheta.*Trig.cphi] +[Aero.X; Aero.Y; Aero.Z]./m+[0;0;Control.Z]./m;

U_dot = v_dot(1);
V_dot = v_dot(2);
W_dot = v_dot(3);


%% P_dot, Q_dot, R_dot
Omega_dot= [ ((In.y-In.z)./In.x).*State.q.*State.r;((In.z-In.x)./In.y).*State.p.*State.r;((In.x-In.y)./In.z).*State.q.*State.q]+[(1/In.x).*Aero.L; (1/In.y).*Aero.M; (1/In.z).*Aero.N]+[(1/In.x).*Control.L; (1/In.y).*Control.M; (1/In.z).*Control.N];

P_dot = Omega_dot(1);
Q_dot = Omega_dot(2);
R_dot = Omega_dot(3);

% %% Ground Collison
% if(State.z>=0)        
%     Z_dot=0; 
%     U_dot=0;
%     V_dot=0;
%     W_dot=0;
%     P_dot=0;
%     Q_dot=0;
%     R_dot=0;
% end

%% Compilation
var_dot=[ X_dot; Y_dot; Z_dot; Phi_dot; Theta_dot; Psi_dot; U_dot; V_dot; W_dot; P_dot; Q_dot; R_dot];


end