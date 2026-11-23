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
tau_grid = linspace(0,1,20); 


K = [4 4];   % vecteur ligne

%% Grille de recherche pour lambda2

lambda_values = linspace(0.8,4,20);
gamma_values = linspace(0.06,0.09901,50);

%% stockage des résultats

best_trace = inf;
best_Q = [];
best_lambda = NaN;
best_gamma = inf;

%% boucle principale

for lambda = lambda_values
  for gamma = gamma_values

     % Matrice fermée
 M = A + alpha1*B*K;
 
 % ----- Variables -----
 Q = sdpvar(2,2);
 
 % ----- Blocs -----
 S = M'*Q + Q*M;
 C = Q*B + lambda*rho1*K';
 E = Q*Bd;
 
 LMI = [ S   C        E;
         C' -2*lambda 0;
         E'  0        -gamma];
 
 % ----- Contraintes -----
 Constraints = [];
 Constraints = [Constraints, Q >= eps_val*eye(2)];
 Constraints = [Constraints, LMI <= -eps_val*eye(4)];
 
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
            best_lambda = lambda;
        end

        if gamma < best_gamma
            best_gamma = gamma;
        end

    end

  end

end

%% résultats

disp('Meilleur résultat trouvé:')
disp(['lambda = ', num2str(best_lambda)])
disp(['gamma = ', num2str(best_gamma)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)