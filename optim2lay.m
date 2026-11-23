clear;
clc;

%% Données du problème

A = [0 1;
     1 0];

B = [0;
    -1];

Bd=[0;
    0.01];

alpha1 = 1/2;
alpha2 = 1;

rho1 = 3/2;
rho2 = 2;

K2 = [4 4];   % vecteur ligne

%% Grilles de recherche

lambda1_values = linspace(0.3,0.45,5);
lambda2_values = linspace(0.2,0.31,5);
tau_grid = linspace(0,1,2); 
gamma_values = linspace(0.095,0.10013,200);

%% stockage des résultats

best_trace = inf;
best_Q = [];
best_lambda1 = NaN;
best_lambda2 = NaN;
best_gamma = inf;

%% boucle principale double

for lambda1 = lambda1_values
 for lambda2 = lambda2_values
    for gamma = gamma_values


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

    M = [M11 M12 M13 Q*Bd;
         M21 M22 M23 0;
         M31 M32 M33 0
         Bd'*Q    0   0  -gamma];

    % contraintes

    Constraints = [];
    Constraints = [Constraints, Q >= 1e-1*eye(2)];
    Constraints = [Constraints, M <= -1e-1*eye(5)];

    % contraintes en tau
    for tau = tau_grid
        a = [tau/3; tau/3];
        LMI_tau = [1 a';
                   a Q];
        Constraints = [Constraints, LMI_tau >= 0];
    end

    % objectif
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
            best_lambda1 = lambda1;
            best_lambda2 = lambda2;
        end

        if gamma < best_gamma
            best_gamma = gamma;
        end

     end

    end
 end
end

%% résultats

disp('Meilleur résultat trouvé:')
disp(['lambda1 = ', num2str(best_lambda1)])
disp(['lambda2 = ', num2str(best_lambda2)])
disp(['gamma = ', num2str(best_gamma)])
disp(['trace(Q) = ', num2str(best_trace)])
disp('Q = ')
disp(best_Q)