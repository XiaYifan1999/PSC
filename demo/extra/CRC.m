function [idx, V, C] = CRC(X, Y, conf)

[N, D] = size(Y);

% Construct kernel matrix K
M = 15;
[U, K] = FourierKernel(X, M);

% Initialization
V = zeros(N,D); 
iter = 1;  
tecr = 1; 
C = zeros(M,D); 
E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); 
gamma = conf.gamma;

%%
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step. 
    E_old = E;
    [P, E] = get_P(Y,V, sigma2 ,gamma, conf.a);

    tecr = abs((E-E_old)/E);

    % M-step. Solve linear system for C.
    P = max(P, conf.minP);
    s = max(sigma2,0.05);
    C = (U'.*repmat(P', [M, 1])*U + conf.lambda*s*K)\(U'.*repmat(P', [M, 1])*Y);

    % Update V and sigma^2
    V = U*C;
    Sp = sum(P);
    sigma2 = sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    % Update gamma
    numcorr = length(find(P > conf.theta));
    gamma = numcorr/size(X, 1);
    
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter = iter + 1;
end

%%
idx = find(P > conf.theta);

%%
function [P, E]=get_P(Y, V, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P = temp1./(temp1+temp2);
E = P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;
