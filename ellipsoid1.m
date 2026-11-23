clc; clear; close all;

Q1 = [1.6 1.4;
      1.4 1.5];

Q2 = [1.01 0.8;
      0.8 0.81];

theta = linspace(0, 2*pi, 2000);
circle = [cos(theta); sin(theta)];

%% Ellipsoïde 1
[V1, D1] = eig(Q1);
A1 = V1 * diag(1./sqrt(diag(D1)));
E1 = A1 * circle;

%% Ellipsoïde 2
[V2, D2] = eig(Q2);
lambda2 = diag(D2);

disp('Valeurs propres de Q2 :')
disp(lambda2)

% Seuil numérique
tol = 1e-10;
lambda2(lambda2 < tol) = tol;

A2 = V2 * diag(1./sqrt(lambda2));
E2 = A2 * circle;

%% Tracé
figure;
plot(E1(1,:), E1(2,:), 'b', 'LineWidth', 2); hold on;
plot(E2(1,:), E2(2,:), 'r', 'LineWidth', 2);

axis equal;
grid on;
%legend('Q_1^{opt}', 'Q_2^{opt}','FontSize',28);
%set(gca,'FontSize',16)
%title('Ellipsoids x^T Q_1^{opt} x = 1 and x^T Q_2^{opt} x = 1 ');



% axes (IMPORTANT : applique aux ticks aussi)
set(gca, ...
    'FontSize',30, ...
    'FontName','Times')

% labels (LaTeX + taille plus grande)
%xlabel('FontSize',26, ...
 %   'FontName','Times')

%ylabel('FontSize',26, ...
 %   'FontName','Times')

% légende (LaTeX + taille cohérente)
hleg = legend({'$Q_1^{opt}$', '$Q_2^{opt}$'}, ...
    'Interpreter','latex', ...
    'FontSize',30, ...
    'Location','best');

% légende (LaTeX + taille cohérente)
%hleg = legend({'continuous NN control','linear control'}, ...
  %  'FontSize',21, ...
  %  'Location','best');

% export EPS (comme dans le script)
set(gcf,'PaperPositionMode','auto')
print('-depsc2','-loose','figure_states.eps')