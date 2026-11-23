clear;
clc;

%% Données du problème

A = [0 1;
     1 0];

B = [0;
    -1];

alpha1 = 1/2;
alpha2 = 1;

rho1 = 3/2;
rho2 = 2;
tau_grid = linspace(0,1,100); 

lambda1 = 0.4;

K2 = [4 4];   % vecteur ligne

%% Grille de recherche pour lambda2

lambda2_values = logspace(-4,1,100);

%% stockage des résultats

best_trace = inf;
best_Q = [];
best_lambda2 = NaN;

%% boucle principale

for lambda2 = lambda2_values

    % variable SDP
    Q = sdpvar(2,2,'symmetric');

    % blocs de la LMI

    M11 = (A + alpha1*alpha2*B*K2)'*Q + Q*(A + alpha1*alpha2*B*K2);

    M12 = alpha2*Q*B + lambda1*rho1*K2';
    M13 = Q*B + lambda2*rho2*alpha1*K2';

    M21 = M12';
    M22 = -2*lambda1;
    M23 = lambda2*rho2;

    M31 = M13';
    M32 = lambda2*rho2;
    M33 = -2*lambda2;

    M = [M11 M12 M13;
         M21 M22 M23;
         M31 M32 M33];

    % contraintes

    Constraints = [];
    Constraints = [Constraints, Q >= 1e-1*eye(2)];
    Constraints = [Constraints, M <= -1e-1*eye(4)];

    % objectif

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
            best_lambda2 = lambda2;
        end

    end

end

%% résultats

disp('Meilleur résultat trouvé:')
disp(['lambda2 = ', num2str(best_lambda2)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)