clc
clear
close all

%% Matrices
A = [0 1;
     1 0];

B = [0;
    -1];

Q = [1.6 1.4;
     1.4 1.5];

%% saturation
sat = @(z) max(-1,min(1,z));

%% contrôle
K=[2,2];

%% dynamique
dynamics = @(t,x) A*x + B*K*x;

%% Condition initiale dans Pi1 et proche du bord
x0 = [1;3.5];

%disp('Verification x0^T Q x0 = ')
%disp(x0'*Q*x0)

%% Simulation
tspan = [0 10];
[t,x] = ode45(dynamics,tspan,x0);

%% Evolution temporelle
figure
plot(t,x(:,1),'LineWidth',2)
hold on
plot(t,x(:,2),'LineWidth',2)

grid on
box on

%% Portrait de phase
%figure
%plot(x(:,1),x(:,2),'LineWidth',2)
%grid on
%box on
%xlabel('x_1')
%ylabel('x_2')
%title('Portrait de phase')

% axes (IMPORTANT : applique aux ticks aussi)
set(gca, ...
    'FontSize',30, ...
    'FontName','Times')

% labels (LaTeX + taille plus grande)
xlabel('$t$', ...
    'Interpreter','latex', ...
    'FontSize',30, ...
    'FontName','Times')

ylabel('$x(t)$', ...
    'Interpreter','latex', ...
    'FontSize',30, ...
    'FontName','Times')

% légende (LaTeX + taille cohérente)
hleg = legend({'$x_1$','$x_2$'}, ...
    'Interpreter','latex', ...
    'FontSize',30, ...
    'Location','best');

% export EPS (comme dans le script)
set(gcf,'PaperPositionMode','auto')
print('-depsc2','-loose','figure_states.eps')