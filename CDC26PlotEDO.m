clc;
clear;
close all;

%% Matrices du système
A = [0 1;
     1 0];

B = [0;
    -1];

Q = [2 1;
     1 2];

%% Fonction saturation
sat = @(z) max(-1, min(1, z));

%% Calcul du contrôle u(x)
function u = control_input(x)
    sat = @(z) max(-1, min(1, z));
    integrand = @(tau) sat((1/3)*tau*(x(1)+x(2)));
    u = 24 * integral(integrand, 0, 1);
end

%% Dynamique du système
function dx = system_dynamics(t, x, A, B)
    u = control_input(x);
    dx = A*x + B*u;
end

%% Génération d'une condition initiale dans Pi1
% x0'Qx0 <= 1
while true
    x0 = 2*rand(2,1) - 1; % tirage aléatoire dans [-1,1]
    if x0'*Q*x0 <= 1
        break;
    end
end

%% Simulation
tspan = [0 10];
[t, x] = ode45(@(t,x) system_dynamics(t,x,A,B), tspan, x0);

%% Affichage des états
figure;
plot(t, x(:,1), 'LineWidth', 2);
hold on;
plot(t, x(:,2), 'LineWidth', 2);
grid on;
xlabel('Times');
ylabel('States');
legend('x_1','x_2');
title('Evolution of the states');

%% Portrait de phase
%figure;
%plot(x(:,1), x(:,2), 'LineWidth', 2);
%grid on;
%xlabel('x_1');
%ylabel('x_2');
%title('Portrait de phase');