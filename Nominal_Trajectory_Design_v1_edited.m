%
% Nominal trajectory computation
%
% Hipotesis 10-15 min fligt; Moon observer is inertial one; Flat Moon, no
% curvature retained; No Moon motion, never rotation retained. Only Moon
% gravity retained; contant gravity considered (take into account the real
% gravity variation in the las 10-15 minutes of the operation), nor
% atmosphere neiter other contact actions on the vehicle. Dynamics on the
% plane; considered cdm trajectory in the plane and attitude dynamics in
% the normal direction to the trajectory plane, that is onlu pitch dynamics.
% 
%
% Reference axes; Cartesian Coordinates, Origin on the Moon surface, x-axis
% along the Moon surface; y-axis along the normal to the Moon surface
% upwards.
%
% Model parameters
%
% g_surface_moon, x_dot_f = 0.0
%
% Design parameters, to be selected, independent variables
%
% x_dot_cero, y_dot_cero, tf_t0, K (pitch attitude parameter tg(theta_N)
%
% Design Constraints, dependent variables to be computed as function of
% independent variables
%
% xf_x0, yf_y0, y_dot_f
%
% Other constraints to be considered
%
% Pitch Attitude dy_dx_1_4_tf_t0, dy_dx_1_2_tf_t0, dy_dx_3_4_tf_t0, dy_dx_7_8_tf_t0
%
clear all
close all
clc
%
rad_to_deg = 180.0/pi;
%
% Parameters
%
g_surface_moon = 1.617; % (m/sec^2)
x_dot_f = 0.0; % (m/sec)
%
% min, delta and numbers of values to be considered for each design
% parameter
%
% x_dot_0 (relative to moon surface)
%
min_x_dot_0 = 1660; % (m/sec)
delta_x_dot_0 = 0.1e00; % (m/sec)
num_x_dot_0 = 10.0e00; % (-)
%
% y_dot_0 (relative to moon surface)
%
min_y_dot_0 = -20.7; % (m/sec)
delta_y_dot_0 = -0.1e00; % (m/sec)
num_y_dot_0 = 11.0e00; % (-)
%
% tf_t0
%
min_tf_t0 = 580; % (sec) 
delta_tf_t0 = 1.0e00 ; % (sec)
num_tf_t0 = 5.0e00; % (-)
%
% Base Thrust Vertical Acceleration
%
min_thrust_vert_acc = 1.6e00; % (m/sec/sec) % 0.9 (m/sec/sec)
delta_thrust_vert_acc = 0.001e00; % (m/sec/sec)
num_thrust_vert_acc = 3.0e00; % (m/sec/sec)
%
%
for index_x_dot_0 = 1:1:num_x_dot_0
%    
    x_dot_0(index_x_dot_0) = min_x_dot_0 + (index_x_dot_0 - 1.0)*delta_x_dot_0;
    xd0 = x_dot_0(index_x_dot_0);
%    
    for index_y_dot_0 = 1:1:num_y_dot_0
%    
        y_dot_0(index_y_dot_0) = min_y_dot_0 + (index_y_dot_0 - 1.0)*delta_y_dot_0;
        yd0 = y_dot_0(index_y_dot_0);
%
        for index_tf_t0 = 1:1:num_tf_t0
%    
            tf_t0(index_tf_t0) = min_tf_t0 + (index_tf_t0 - 1.0)*delta_tf_t0;
            tf0 = tf_t0(index_tf_t0);
%
            for index_thrust_vert_acc = 1:1:num_thrust_vert_acc
%    
                thrust_vert_acc(index_thrust_vert_acc) = min_thrust_vert_acc + (index_thrust_vert_acc - 1.0)*delta_thrust_vert_acc;
                thrust_vert_acc_c = thrust_vert_acc(index_thrust_vert_acc);
%
% xf_x0
                xf_x0(index_x_dot_0,index_tf_t0,index_y_dot_0,index_thrust_vert_acc) = xd0*tf0 + 0.5e00*(x_dot_f - xd0)/tf0*tf0*tf0; % (m) x speed down
%
% yf_y0
                yf_y0(index_y_dot_0,index_thrust_vert_acc,index_x_dot_0,index_tf_t0) = yd0*tf0 + 0.5e00*(-g_surface_moon + thrust_vert_acc_c)*tf0*tf0; % (m) y speed up
%
% y_dot_f
                y_dot_f(index_y_dot_0,index_thrust_vert_acc,index_x_dot_0,index_tf_t0) = yd0 + (-g_surface_moon + thrust_vert_acc_c)*tf0; % (m/sec) y speed up
%
% Pitch Attitude Parameter

                par_att = sqrt((x_dot_f - xd0)/tf0*(x_dot_f - xd0)/tf0 + thrust_vert_acc_c*thrust_vert_acc_c);
                theta_pitch_att(index_x_dot_0,index_y_dot_0,index_tf_t0,index_thrust_vert_acc) = rad_to_deg*atan2(thrust_vert_acc_c/par_att,(x_dot_f - xd0)/tf0/par_att); % (deg
                                
% dy_dx_1_4_tf_t0
%
                par_t = 0.25e00*tf0; % (sec)
                d_x = xd0 + (x_dot_f - xd0)/tf0*par_t; % (m) x speed down
                d_y = yd0 + (-g_surface_moon + thrust_vert_acc_c)*par_t; % (m) y speed up
                dy_dx_1_4_tf_t0(index_x_dot_0,index_y_dot_0,index_tf_t0,index_thrust_vert_acc) = d_y/d_x; % (-)
%
% dy_dx_1_2_tf_t0
%
                par_t = 0.5e00*tf0; % (sec)
                d_x = xd0 + (x_dot_f - xd0)/tf0*par_t; % (m) x speed down
                d_y = yd0 + (-g_surface_moon + thrust_vert_acc_c)*par_t; % (m) y speed up
                dy_dx_1_2_tf_t0(index_x_dot_0,index_y_dot_0,index_tf_t0,index_thrust_vert_acc) = d_y/d_x; % (-)
%
% dy_dx_3_4_tf_t0
%
                par_t = 0.75e00*tf0; % (sec)
                d_x = xd0 + (x_dot_f - xd0)/tf0*par_t; % (m) x speed down
                d_y = yd0 + (-g_surface_moon + thrust_vert_acc_c)*par_t; % (m) y speed up
                dy_dx_3_4_tf_t0(index_x_dot_0,index_y_dot_0,index_tf_t0,index_thrust_vert_acc) = d_y/d_x; % (-)
%
% dy_dx_7_8_tf_t0
%
                par_t = (7.0e00/8.0e00)*tf0; % (sec)
                d_x = xd0 + (x_dot_f - xd0)/tf0*par_t; % (m) x speed down
                d_y = yd0 + (-g_surface_moon + thrust_vert_acc_c)*par_t; % (m) y speed up
                dy_dx_7_8_tf_t0(index_x_dot_0,index_y_dot_0,index_tf_t0,index_thrust_vert_acc) = d_y/d_x; % (-)
%
            end
%
        end
%
    end
%
end
%
%% Graphical output, every table s function of x_dot_0 and y_dot_0 with
% parameters tf_t0 and K
%
i_fig = 0;

%% xf_x0

i_fig = i_fig + 1;

label_title = ['Desplazamiento en x (%) vs velocidad X inicial y  tiempo de maniobra total. No tiene parametros porque los canales estan desacoplados.'];
label_z = ['Desplazamiento en X (%)'];
label_y = ['Velocidad inicial en X (m/s)'];
label_x = ['Tiempo de maniobra total (s)'];

% Da igual que valor de velocidad incial y ponemos y de aceleración vertical
% ya que en este modelo, son independientes los canales.

figure(i_fig)
surf(tf_t0,x_dot_0,(xf_x0(:,:,1,1)-481520)/481520*100,'FaceAlpha',.3,'EdgeAlpha',.3)
grid on
title(label_title)
zlabel(label_z)
xlabel(label_x)
ylabel(label_y)

% Para este calculo, los resultados más apropiados son: t=580s y
% x_dot_0=1660m/s

%% yf_y0

i_fig = i_fig + 1;

label_title = ['Desplazamiento en y (%) vs velocidad Y inicial y  ac. vertical cte. Tiempo=580s'];
label_z = ['Desplazamiento en Y (%)'];
label_y = ['Velocidad inicial en Y (m/s)'];
label_x = ['Aceleración vertical de empuje (m/s^2)'];

figure(i_fig)
surf(thrust_vert_acc,y_dot_0,(yf_y0(:,:,1,1)+15084)/15084*100,'FaceAlpha',.3,'EdgeAlpha',.3)
grid on
title(label_title)
zlabel(label_z)
xlabel(label_x)
ylabel(label_y)

% Para este calculo, los resultados más apropiados son: y_dot_0=-21.1m/s y
% vert_thrust=1.6m/s^2

%% y_dot_f

%y_dot_final = -21.1 + (-g_surface_moon + 1.6)*580

i_fig = i_fig + 1;

label_title = ['Velocidad final en Y (m/s) vs velocidad Y inicial y  ac. vertical cte. Tiempo=580s'];
label_z = ['Velocidad final en Y (m/s)'];
label_y = ['Velocidad inicial en Y (m/s)'];
label_x = ['Aceleración vertical de empuje (m/s^2)'];

figure(i_fig)
surf(thrust_vert_acc,y_dot_0,y_dot_f(:,:,1,1),'FaceAlpha',.3,'EdgeAlpha',.3)
grid on
title(label_title)
zlabel(label_z)
xlabel(label_x)
ylabel(label_y)





