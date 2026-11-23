clear;
clc;
format long
digits(8)

%% Données du problème

A = [0 1;
     1 0];

B = [0;
    -1];

alpha1 = 1/2;
alpha2 = 1;

rho1 = 3/2;
rho2 = 2;
eps_val = 1e-1;
tau_grid = linspace(0,1,100); 


K = [4 4];   % vecteur ligne

%% Grille de recherche pour lambda2

lambda_values = linspace(0.5,1,50);

%% stockage des résultats

best_trace = inf;
best_Q = [];
best_lambda = NaN;

%% boucle principale

for lambda = lambda_values

     % Matrice fermée
 M = A + alpha1*B*K;
 
 % ----- Variables -----
 Q = sdpvar(2,2);
 
 % ----- Blocs -----
 S = M'*Q + Q*M;
 C = Q*B + lambda*rho1*K';
 
 LMI = [ S C;
 C' -2*lambda ];
 
 % ----- Contraintes -----
 Constraints = [];
 Constraints = [Constraints, Q >= eps_val*eye(2)];
 Constraints = [Constraints, LMI <= -eps_val*eye(3)];
 
 % borne pour éviter problème non borné
  Constraints = [Constraints, trace(Q) >=1e-6];
 
 % ----- Objectif : minimiser trace(Q) -----
 Objective = trace(Q);
    % résolution

    options = sdpsettings('solver','sdpt3','verbose',0);

    for tau = tau_grid
 
     a = [tau/3; tau/3]; % vecteur colonne
 
     LMI_tau = [1 a';
                  a Q];
 
     Constraints = [Constraints, LMI_tau >= 0];
 end

    sol = optimize(Constraints,Objective,options);

    if sol.problem == 0

        Q_val = value(Q);
        tr = trace(Q_val);

        if tr < best_trace
            best_trace = tr;
            best_Q = Q_val;
            best_lambda = lambda;
        end

    end

end

%% résultats

disp('Meilleur résultat trouvé:')
disp(['lambda = ', num2str(best_lambda)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)