% Nominal Landing trajectory for the selected design parameters version v1

% Horizontal acceleration selected to verify the final point hovering condition
% Thrust vertical acceleration selected as input parameter
% Pitch Attitude computed

close all
clear all
clc

global tf_t0 thrust_vert_acc x_dot_0 y_dot_0 x_dot_f g_surface_moon
global delta_g delta_x_0 delta_x_dot_0 delta_y_0 delta_y_dot_0 delta_thrust_vert_acc sigma_noise_x sigma_noise_y v_time_noise noise_x noise_y
%
rad_to_deg = 180.0/pi;

% Input data

% tf_t0 (sec), thrust_vert_acc (m/sec/sec), x_dot_0 (m/sec), y_dot_0 (m/sec), x_dot_f = 0 (m/sec)

tf_t0 = 580; % (sec)
thrust_vert_acc = 1.6; % (m/sec/sec)
x_dot_0 = 1660; % (m/sec)
y_dot_0 = -21.1; % (m/sec)
x_dot_f = 0.0; % (m/sec)

% Model parameter

g_surface_moon = 1.617; % (m/sec^2)

% Time step

Time_step = 1.0e-02; % (sec)

% Time vector

v_t = [0:Time_step:tf_t0];

% Nominal Trajectory

%% Analitical solution

x_x0 = x_dot_0*v_t + 0.5e00*(x_dot_f - x_dot_0)/tf_t0*v_t.*v_t;
x_dot = x_dot_0 + (x_dot_f - x_dot_0)/tf_t0*v_t;
y_y0 = y_dot_0*v_t + 0.5e00*(-g_surface_moon + thrust_vert_acc)*v_t.*v_t;
y_dot = y_dot_0 + (-g_surface_moon + thrust_vert_acc)*v_t;
provi = sqrt((x_dot_f - x_dot_0)/tf_t0*(x_dot_f - x_dot_0)/tf_t0 + thrust_vert_acc*thrust_vert_acc); 
pitch_angle_ana = rad_to_deg*atan2(thrust_vert_acc/provi,(x_dot_f - x_dot_0)/tf_t0/provi);

for index_barrido = 1:1:max(size(x_dot))
    provi = sqrt(x_dot(index_barrido)*x_dot(index_barrido) + y_dot(index_barrido)*y_dot(index_barrido));
    flight_path_angle_ana(index_barrido) = rad_to_deg*atan2(y_dot(index_barrido)/provi,x_dot(index_barrido)/provi);
end


%% Numerical solution (using ode45, Runge-Kutta)

% Initial condition

trajectory_num_ini = [0.0 x_dot_0 0.0 y_dot_0];

% Dynamic evolution

[t,trajectory_num] = ode45('der_trajectory_v1',[0 tf_t0],trajectory_num_ini);
provi = sqrt((x_dot_f - x_dot_0)/tf_t0*(x_dot_f - x_dot_0)/tf_t0 + thrust_vert_acc*thrust_vert_acc); 
pitch_angle_num = rad_to_deg*atan2(thrust_vert_acc/provi,(x_dot_f - x_dot_0)/tf_t0/provi);

for index_barrido = 1:1:max(size(trajectory_num))
    provi = sqrt(trajectory_num(index_barrido,2)*trajectory_num(index_barrido,2) + trajectory_num(index_barrido,4)*trajectory_num(index_barrido,4));
    flight_path_angle_num(index_barrido) = rad_to_deg*atan2(trajectory_num(index_barrido,4)/provi,trajectory_num(index_barrido,2)/provi);
end

%% Perturbed Trajectory

% Perturbation Parameters

delta_g = 1.0e-01; % (m/sec/sec)
delta_x_0 = 1.0e03; % (m)
delta_x_dot_0 = 1.0e00; % (m/sec)
delta_y_0 = -1.0e02; % (m)
delta_y_dot_0 = -1.0e00; % (m/sec)
delta_thrust_vert_acc = 6.0e-02; % (m/sec)

% Perturbation percentages

% delta_g/g_surface_moon*100
% delta_x_0/481520*100
% delta_x_dot_0/x_dot_0*100
% delta_y_0/-15084*100
% delta_y_dot_0/y_dot_0*100
% delta_thrust_vert_acc/thrust_vert_acc*100

% Noise effect

sigma_noise_x = 1.0e-01; % (m/sec/sec)
sigma_noise_y = 1.0e-01; % (m/sec/sec)
v_time_noise = [0:Time_step:tf_t0];

rng(7);

noise_x = sigma_noise_x*randn(size(v_time_noise));
noise_y = sigma_noise_y*randn(size(v_time_noise));

% Numerical solution (using ode45, Runge-Kutta)

% Initial condition

trajectory_num_pert_ini = [delta_x_0 x_dot_0+delta_x_dot_0 delta_y_0 y_dot_0+delta_y_dot_0];

% Dynamic evolution

[t_pert,trajectory_num_pert] = ode45('der_trajectory_pert_v1',[0 tf_t0],trajectory_num_pert_ini);
provi = sqrt((x_dot_f - x_dot_0)/tf_t0*(x_dot_f - x_dot_0)/tf_t0 + thrust_vert_acc*thrust_vert_acc); 
pitch_angle_num_pert = rad_to_deg*atan2(thrust_vert_acc/provi,(x_dot_f - x_dot_0)/tf_t0/provi);

for index_barrido = 1:1:max(size(trajectory_num_pert))
    provi = sqrt(trajectory_num_pert(index_barrido,2)*trajectory_num_pert(index_barrido,2) + trajectory_num_pert(index_barrido,4)*trajectory_num_pert(index_barrido,4));
    flight_path_angle_num_pert(index_barrido) = rad_to_deg*atan2(trajectory_num_pert(index_barrido,4)/provi,trajectory_num_pert(index_barrido,2)/provi);
end

%% Graphical Output

i_fig = 0.0;
i_fig = i_fig + 1;

s10 = ['Apollo Moon Planar Landing Nominal Trayectories; Flight Time = ' num2str(tf_t0) ' (sec); Initial x-velocity = ' num2str(x_dot_0) ' (m/sec); Initial y-velocity = ' num2str(y_dot_0) ' (m/sec)'];

figure(i_fig)

% X Position

subplot(211)
plot(t,trajectory_num(:,1),'g','Linewidth',1.5)
hold on
plot(t_pert,trajectory_num_pert(:,1),'r','Linewidth',1.5)
hold on
plot(v_t,x_x0,'b','Linewidth',1.5)
axis([min(t) max(t) min(trajectory_num(:,1)) max(trajectory_num(:,1))])
grid
legend('Posición X numérica','Posición X analitica', 'Posición X perturbada')
title(s10)
ylabel('X Moon Reference Frame (m)')

% X Speed

subplot(212)
plot(t,trajectory_num(:,2),'g','Linewidth',1.5)
hold on
plot(t_pert,trajectory_num_pert(:,2),'r','Linewidth',1.5)
hold on
plot(v_t,x_dot,'b','Linewidth',1.5)
axis([min(t) max(t) min(trajectory_num(:,2)) max(trajectory_num(:,2))])
grid
legend('Velocidad X numérica','Velocidad X analitica', 'Velocidad X perturbada')
xlabel('Time (sec)')
ylabel('X velocity Moon Reference Frame (m/sec)')

i_fig = i_fig + 1;

figure(i_fig)

% Y Position

subplot(211)
plot(t,trajectory_num(:,3),'g','Linewidth',1.5)
hold on
plot(v_t,y_y0,'r','Linewidth',1.5)
hold on
plot(t_pert,trajectory_num_pert(:,3),'b','Linewidth',1.5)
axis([min(t) max(t) min(trajectory_num(:,3)) max(trajectory_num(:,3))])
grid
legend('Posición Y numérica','Posición Y analitica', 'Posición Y perturbada')
title(s10)
ylabel('Y Moon Reference Frame (m)')

% Y Speed

subplot(212)
plot(t,trajectory_num(:,4),'g','Linewidth',1.5)
hold on
plot(v_t,y_dot,'r','Linewidth',1.5)
hold on
plot(t_pert,trajectory_num_pert(:,3),'b','Linewidth',1.5)
axis([min(t) max(t) min(trajectory_num(:,4)) max(trajectory_num(:,4))])
grid
legend('Velocidad Y numérica','Velocidad Y analitica', 'Velocidad Y perturbada')
xlabel('Time (sec)')
ylabel('Y velocity Moon Reference Frame (m/sec)')

% Trajectories

i_fig = i_fig + 1;

figure(i_fig)
plot(trajectory_num(:,1),trajectory_num(:,3),'g','Linewidth',1.5)
hold on
plot(x_x0,y_y0,'r','Linewidth',1.5)
hold on
plot(trajectory_num_pert(:,1),trajectory_num_pert(:,3),'b','Linewidth',1.5)
axis([min(trajectory_num(:,1)) max(trajectory_num(:,1)) min(trajectory_num(:,3)) max(trajectory_num(:,3))])
grid
legend('Trayectoria numérica','Trayectoria analitica', 'Trayectoria perturbada')
title(s10)
xlabel('X Moon Reference Frame (m)')
ylabel('Y Moon Reference Frame (m)')

i_fig = i_fig + 1;

figure(i_fig)

% Pitch Angle

subplot(211)
plot(t,pitch_angle_num*ones(1,max(size(t))),'g','Linewidth',1.5)
hold on
plot(v_t,pitch_angle_ana*ones(1,max(size(v_t))),'r','Linewidth',1.5)
hold on
plot(t_pert,pitch_angle_num_pert*ones(1,max(size(t_pert))),'b','Linewidth',1.5)
grid
title(s10)
ylabel('Pitch angle (º)')
legend('Numerical','Analitical','Perturbed')

% Flight Path Angle

subplot(212)
plot(t,flight_path_angle_num,'g','Linewidth',1.5)
hold on
plot(v_t,flight_path_angle_ana,'r','Linewidth',1.5)
hold on
plot(t_pert,flight_path_angle_num_pert,'b','Linewidth',1.5)
grid
legend('Numerical','Analitical','Perturbed')
xlabel('Time (sec)')
ylabel('Flight Path angle (º)')

%% Almacenamos la trayectoria nominal para alimentar el control

t_nom = t;
trajectory_num_nom = trajectory_num;

save apollo_moon_landing_v1_case_1 t_nom trajectory_num_nom pitch_angle_num flight_path_angle_num
