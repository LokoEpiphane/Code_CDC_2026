clear;
clc;
format long
digits(8)

%% Données du problème

A = [0 1;
     1 0];

B = [0;
    -1];

Bd =[0;
     0.01];

alpha1 = 1/2;
alpha2 = 1;

rho1 = 3/2;
rho2 = 2;
eps_val = 1e-1;


K = [4 4];   % vecteur ligne

%% Grille de recherche pour lambda2

gamma_values = linspace(0.07,0.12,3500);

%% stockage des résultats

best_trace = inf;
best_Q = [];
best_gamma = inf;

%% boucle principale

  for gamma = gamma_values

     % Matrice fermée
 M = A + alpha1*B*K;
 
 % ----- Variables -----
 Q = sdpvar(2,2);
 
 % ----- Blocs -----
 S = M'*Q + Q*M;
 E = Q*Bd;
 
 LMI = [ S   E;
         E'  -gamma];
 
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

    sol = optimize(Constraints,Objective,options);

    if sol.problem == 0

        Q_val = value(Q);
        tr = trace(Q_val);

        if tr < best_trace
            best_trace = tr;
            best_Q = Q_val;
        end

        if gamma < best_gamma
            best_gamma = gamma;
        end

    end

  end

%% résultats

disp('Meilleur résultat trouvé:')
disp(['gamma = ', num2str(best_gamma)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)