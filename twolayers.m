clear; clc;

%% Données du problème
A = [0 1;
     1 0];

B = [0;
     -1];

alpha1 = 1/2;
alpha2 = 1/2;
rho1   = 3/2;
rho2   = 3/2;
lambda1 = 1;

e = [1;1];

%% Grille de valeurs non nulles de w3
w3_values = linspace(60,200,200);
w3_values(w3_values==0) = [];
tau_grid = linspace(0,1,20); 

best_trace = inf;
best_Q = [];
best_lambda2 = [];
best_w3 = [];

for w3 = w3_values
    
    K = (w3/18)*[1 1];
    
    % Variables de décision
    Q = sdpvar(2,2);
    lambda2 = sdpvar(1);
    
    % Terme dynamique
    Acl = A + alpha1*alpha2*B*K;
    
    % Construction des blocs
    M11 = Acl'*Q + Q*Acl;
    
    M12 = (alpha2*w3/2)*Q*B + (lambda2*rho1/6)*e;
    M13 = w3*Q*B + (alpha1*lambda2*rho2/12)*e;
    
    M22 = -2*lambda2;
    M23 = (lambda2*rho2)/2;
    M33 = -2*lambda2;
    
    % Matrice LMI complète
    M = [ M11     M12        M13;
          M12'    M22        M23;
          M13'    M23        M33 ];
      
    % Contraintes
    Constraints = [];
    Constraints = [Constraints, Q >= 1e-6*eye(2)];
    Constraints = [Constraints, lambda2 >= 1e-6];
    Constraints = [Constraints, M <= -1e-6*eye(4)];
    
    % Objectif : minimiser trace(Q)
    Objective = trace(Q);
    
    options = sdpsettings('solver','sdpt3','verbose',0);

    %for tau = tau_grid
 
        %a = [tau/3; tau/3]; % vecteur colonne
 
        %LMI_tau = [1 a';
        %             a Q];
 
        %Constraints = [Constraints, LMI_tau >= 0];
    %end
    diagnostics = optimize(Constraints,Objective,options);
    
    if diagnostics.problem == 0
        Qval = value(Q);
        lambda2val = value(lambda2);
        tr = trace(Qval);
        
        if tr < best_trace
            best_trace = tr;
            best_Q = Qval;
            best_lambda2 = lambda2val;
            best_w3 = w3;
        end
    end
end

%% Résultat
if isempty(best_Q)
    disp('Aucune solution trouvée.')
else
    disp('===== Meilleure solution trouvée =====')
    disp(['w3 optimal = ', num2str(best_w3)])
    disp('Q = ')
    disp(best_Q)
    disp(['lambda2 = ', num2str(best_lambda2)])
    disp(['trace(Q) = ', num2str(best_trace)])
end
