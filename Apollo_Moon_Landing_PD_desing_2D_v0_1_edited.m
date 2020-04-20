clear all
close all
clc

%% INICIALIZACION DE VARIABLES

% Matriz Dynamic state vector [x x_dot y y_dot]

A_free = [0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0;]; % System Matrix

% Free Eigen structure

[V,D] = eig(A_free);
[W,J] = jordan(A_free);

% Genero el rango de pulsaciones (frecuencias) que voy a estudiar espaciando puntos logaritmicamente entre las decadas 10^d1 y 10^d2 (rad/sec)

N = 500.;
d_1 = -2.0;
d_2 = 2.0;

v_w = logspace(d_1,d_2,N);

% Definicion de amortiguamiento y pulsacion natural del sistema original

double_epsilon_omega_n_x = 0.0; % El sistema natural no tiene amortiguamiento
omega_n_square_x = 0.0; % El sistema natural no tiene pulsacion natural

% Definicion de la playa de ganancias que se van a analizar

K_p_x_ini = 0.005; % La minima consigue polos nulos
K_p_x_fin = 0.02; % 2.5e-00
num_K_p_x = 4; % 1.5e+01
delta_K_p_x = (K_p_x_fin - K_p_x_ini)/num_K_p_x;
%
K_d_x_ini = 0.15; % la minima sin amortiguamiento adicional
K_d_x_fin = 0.5; % 
num_K_d_x = 4;
delta_K_d_x = (K_d_x_fin - K_d_x_ini)/num_K_d_x; 


%% CALCULOS

for i = 1:1:num_K_d_x+1 % Bucle en K_d

    % Vector de Kds
    
    K_d_x = K_d_x_ini + (i - 1.0)*delta_K_d_x;
    v_K_d_x(i) = K_d_x;
    
    % Caracteristicas sistema controlado
  
    double_epsilon_omega_c_x(i) = double_epsilon_omega_n_x + K_d_x;
    
    for j = 1:1:num_K_p_x+1 % Bucle en K_p
        
        % Vector de Kps
        
        K_p_x = K_p_x_ini + (j - 1.0)*delta_K_p_x;
        v_K_p_x(j) = K_p_x;
        labels(j)="Kp="+K_p_x;
    
        % Caracteristicas sistema controlado
  
        omega_c_square_x(j)= omega_n_square_x + K_p_x;

        % Sistema en Bucle Abierto
        %%% Funcion de transferencia del sistema G(s)=Kp/(s*s + Kd*s)
    
        ft_x(i,j) = zpk([],[0.,-K_d_x],[K_p_x]);
        
        % Sistema en bucle cerrado
        
        sys(i,j)=feedback(ft_x(i,j),1);

        % Calculo los margenes de fase y de ganancia, y las correspondientes frecuencia de margenes de fase y de ganancia

        [margen_ganancia_x(i,j),margen_fase_x(i,j),frec_ganancia_x(i,j),frec_fase_x(i,j)] = margin(ft_x(i,j));
        
        omega_c_x(j) = sqrt(omega_n_square_x + K_p_x); 
        epsilon_c_x(i,j) = double_epsilon_omega_c_x(i)/2.0/omega_c_x(j);
        
        % Error escalon
        
        error_escalon_x(j) = K_p_x/(omega_c_x(j)^2) - 1.0 ;
        error_escalon_x_dot(i,j) = K_d_x/(omega_c_x(j)^2) - 1.0 ;
        
        % SOBRE O SUBAMORTIGUADO
        
        if epsilon_c_x(i,j)<1.0
            
            caso(i,j)=1.;
            
            % Pulsacion forzada
            
            omega_c_d_x(i,j) = omega_c_x(j)*sqrt(1.0-epsilon_c_x(i,j)*epsilon_c_x(i,j));
            
            % Tiempo de asentamiento a un k por ciento del valor estacionario tras el escalon
            
            k = 5.0;
            t_a_k_x(i,j) = -log((k/100.)*sqrt(1.0-epsilon_c_x(i,j)*epsilon_c_x(i,j)))/epsilon_c_x(i,j)/omega_c_x(j);
        else
            caso(i,j)=0.;
            omega_c_d_x(i,j) = 0;
            t_a_k_x(i,j)=0;
        end
        
    end 
    
end

%% GRAFICOS

n_g_x = 0.;

% DIAGRAMAS DE BODE

% n_g_x = n_g_x + 1.;
% figure(n_g_x)
% title(['Diagrama de Bode de planta en open loop with selected K d x = ' num2str(v_K_d_x(i)) ' X channel']) %  
% grid

for i = 1:1:num_K_d_x+1
    
    n_g_x = n_g_x + 1.;
    figure(n_g_x)
    
    for j = 1:1:num_K_p_x+1
    bode(ft_x(i,j),v_w)
    hold on
    end
    
    title(['Diagrama de Bode de planta en open loop with selected K d x = ' num2str(v_K_d_x(i)) ' X channel']) 
    legend(labels)
    grid
    
end 

for i = 1:1:num_K_d_x+1
    
    n_g_x = n_g_x + 1.;
    figure(n_g_x)
    
    for j = 1:1:num_K_p_x+1
    bode(sys(i,j),v_w)
    hold on
    end
    
    title(['Diagrama de Bode de planta en closed loop with selected K d x = ' num2str(v_K_d_x(i)) ' X channel']) 
    legend(labels)
    grid
    
end 

% POLOS DE LOS SISTEMAS

% for i = 1:1:num_K_d_x+1
%     
%     n_g_x = n_g_x + 1.;
%     figure(n_g_x)
%     
%     for j = 1:1:num_K_p_x+1
%     pzmap(sys(i,j))
%     hold on
%     end
%     
%     title(['Polos con Kd = ' num2str(v_K_d_x(i)) ' X channel']) 
%     legend(labels)
%     grid
%     
% end 

% DIAGRAMAS DE NYQUIST

% n_g_x = n_g_x + 1.;
% figure(n_g_x)
% title(['Diagrama de Nyquist de planta en open loop with selected K d x = ' num2str(v_K_d_x(i)) ' X channel']) %  
% grid

% for i = 1:1:num_K_d_x+1
%     
%     n_g_x = n_g_x + 1.;
%     figure(n_g_x)
%     
%     for j = 1:1:num_K_p_x+1
%     nyquist(ft_x(i,j),v_w)
%     hold on
%     end
%     
%     title(['Diagrama de Nyquist de planta en open loop with selected K d x = ' num2str(v_K_d_x(i)) ' X channel']) 
%     legend(labels)
%     grid
%     
% end 

% MARGEN DE GANANCIA Y MARGEN DE FASE CON SUS RESPECTIVAS FASES

% n_g_x = n_g_x + 1.;
% figure(n_g_x)
% 
% label_title = ['Margen de ganancia respecto a Kp y Kd'];
% label_z = ['Margen de Ganancia (dB)'];
% label_y = ['Kd'];
% label_x = ['Kp'];
% 
% surf(v_K_p_x,v_K_d_x,margen_ganancia_x,'FaceAlpha',.3,'EdgeAlpha',.3)
% grid on
% title(label_title)
% zlabel(label_z)
% xlabel(label_x)
% ylabel(label_y)

% n_g_x = n_g_x + 1.;
% figure(n_g_x)
% 
% label_title = ['Margen de fase respecto a Kp y Kd'];
% label_z = ['Margen de Fase (º)'];
% label_y = ['Kd'];
% label_x = ['Kp'];
% 
% surf(v_K_p_x,v_K_d_x,margen_fase_x,'FaceAlpha',.3,'EdgeAlpha',.3)
% grid on
% title(label_title)
% zlabel(label_z)
% xlabel(label_x)
% ylabel(label_y)

% RESPUESTA A ESCALON 

for i = 1:1:num_K_d_x+1
    
    n_g_x = n_g_x + 1.;
    figure(n_g_x)
    
    for j = 1:1:num_K_p_x+1
    step(sys(i,j))
    hold on
    end
    
    title(['Respuestas a impulso con K d x = ' num2str(v_K_d_x(i)) ' X channel']) 
    legend(labels)
    grid
    
end 

% ERROR A ESCALON X Y ERROR A ESCALON X DOT

n_g_x = n_g_x + 1.;
figure(n_g_x)

%subplot(211)
plot(v_K_p_x,error_escalon_x,'LineWidth',1.3)
title(['Error escalon X segun Kp']) 
xlabel('Kp')
ylabel('Error escalon X')
grid on

n_g_x = n_g_x + 1.;
figure(n_g_x)

%subplot(212)
surf(v_K_p_x,v_K_d_x,error_escalon_x_dot,'FaceAlpha',.3,'EdgeAlpha',.3)
title(['Error escalon X_dot segun Kp y Kd'])
xlabel('Kp')
ylabel('Kd')
zlabel('Error escalon Velocidad X')
grid on

% PULSACIÓN FORZADA

n_g_x = n_g_x + 1.;
figure(n_g_x)

surf(v_K_p_x,v_K_d_x,omega_c_d_x,'FaceAlpha',.3,'EdgeAlpha',.3)
title(['Pulsacion forzada segun Kp y Kd'])
xlabel('Kp')
ylabel('Kd')
zlabel('Pulsacion (1/s)')
grid on

% TIEMPO DE ASENTAMIENTO A UN 5% 

n_g_x = n_g_x + 1.;
figure(n_g_x)

surf(v_K_p_x,v_K_d_x,t_a_k_x,'FaceAlpha',.3,'EdgeAlpha',.3)
title(['Tiempo de asentamiento a 5% segun Kp y Kd'])
xlabel('Kp')
ylabel('Kd')
zlabel('Tiempo de asentamiento (s)')
grid on

% PLANO KP Y KD Y CASUÍSTICA

% n_g_x = n_g_x + 1.;
% figure(n_g_x)
% 
% surf(v_K_p_x,v_K_d_x,caso,'FaceAlpha',.3,'EdgeAlpha',.3)
% title(['Casuística'])
% xlabel('Kp')
% ylabel('Kd')
% zlabel('-')
% grid on
