clc
clear
close all


%% Matrices
A = [0 1;
     1 0];

B = [0;
    -1];

%% saturation
sat = @(z) max(-1,min(1,z));
K=[0.01,0.1]';

%% contrôle
control_input = @(x) 72.*integral(@(tau) sat(2.*tau.*integral(@(s) s.*sat((1/3).*s.*(x(1)+x(2))),0,1)),0,1);

%% dynamique
dynamics = @(t,x) A*x + B*control_input(x);

%% Condition initiale dans Pi1 et proche du bord
x0 = [1;3.5];


%% Simulation
tspan = [0 10];
[t,x] = ode45(dynamics,tspan,x0);

%% Calcul du contrôle le long du temps
u = zeros(length(t),1);
u_lin = zeros(length(t),1);

for i = 1:length(t)
    xi = x(i,:)';
    u(i) = control_input(xi);
    u_lin(i) = [2 2]*xi;
end

%% Evolution temporelle
figure
plot(t,x(:,1),'LineWidth',2)
hold on
plot(t,x(:,2),'LineWidth',2)

grid on
box on

%% Evolution des contrôles
%figure
%plot(t,u,'LineWidth',2)
%hold on
%plot(t,u_lin,'--','LineWidth',2)

%grid on
%box on



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

% légende (LaTeX + taille cohérente)
%hleg = legend({'continuous NN control','linear control'}, ...
 %   'FontSize',30, ...
 %   'Location','best');

% export EPS (comme dans le script)
set(gcf,'PaperPositionMode','auto')
print('-depsc2','-loose','figure_states.eps')