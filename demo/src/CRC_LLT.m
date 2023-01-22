function [idx, V, C] = CRC_LLT(X, Y, conf, P_init,	Ci)

[N, D] = size(Y);
% Construct kernel matrix K
M = 15;
[U, L] = FourierKernel(X, M);

% Initialization
if isempty(Ci)
    V = zeros(N,D); %
else
    V = U*Ci;
end

iter = 1;  
tecr = 1; 
C = zeros(M,D); 
E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); 
gamma = conf.gamma;
Kn = conf.Kn;
if(Kn>D) 
  tol=1e-5; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end

if Kn>N-1
    Kn=N-1;
end

% tic;

% Compute IsubW
kdtreeX = vl_kdtreebuild(X');
[neighborX, ~] = vl_kdtreequery(kdtreeX, X', X', 'NumNeighbors', Kn+1);
neighborhood = neighborX(2:(1+Kn),:);
W = zeros(Kn,N);%wij represprent xi  is the set of neighbors of xj and the wight is wij
% WW = sparse(1:N,1:N,zeros(1,N),N,N,4*Kn*N);%N*N  
IsubW  = sparse(1:N,1:N,ones(1,N),N,N,4*Kn*N);%N*N
for i = 1:N
    z = X(neighborhood(:,i),:) - repmat(X(i,:),Kn,1); % shift ith pt to origin  K*D
    G = z*z';                                         % local covariance  K*K
    G = G + eye(Kn,Kn)* tol * trace(G);                 % regularlization
    W(:,i) = G\ones(Kn,1);                             % solve Gw = 1
    W(:,i) = W(:,i)/sum(W(:,i));                     % normalize
    w = W(:,i);
    j = neighborhood(:,i);
    % WW(i,j) = w';
    IsubW(i,j) = IsubW(i,j) - w';
end
%%

while (iter < conf.MaxIter) && (tecr > conf.ecr)  && (sigma2 > 1e-8) % 
    % E-step. 
    E_old = E;
    if iter>1 || isempty(P_init)
        [P, E] = get_P(Y,V, sigma2 ,gamma, conf.a);
    else
        P = P_init;[~, E] = get_P(Y,V, sigma2 ,gamma, conf.a);
    end
    Psm = sparse(1:N,1:N,P,N,N,N);
    E = E + conf.lambda1*trace(C'*L*C) + conf.lambda2*norm(sqrt(Psm) * IsubW * V,'fro');
%     disp(['L constraint:',num2str(conf.lambda1*trace(C'*L*C))]);
%     disp(['LLE constraint:', num2str(conf.lambda2/2*norm(sqrt(Psm) * IsubW * V,'fro'))]);
    Q = IsubW'*Psm*IsubW;
    
    tecr = abs((E-E_old)/E);

    % M-step. Solve linear system for C.
    P = max(P, conf.minP);
    s = max(sigma2,0.05);
    C = (U'.*repmat(P', [M, 1])*U + conf.lambda1*s*L+conf.lambda2*s*U'*Q*U)\(U'.*repmat(P', [M, 1])*Y);

    % Update V and sigma^2
    V = U*C;
    Sp = sum(P);
    sigma2 = sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    % Update gamma
    numcorr = length(find(P > conf.tau));
    gamma = numcorr/size(X, 1);
    
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter = iter + 1;
end

%%
idx = find(P > conf.tau);

%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y, V, sigma2 ,gamma, a) 
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-V).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P = temp1./(temp1+temp2);
E = P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;  % 
