function [K, L] = FourierKernel(Xn,c)

N = size(Xn,1);
% ascending order of eigenvalue P_i^2 + Q_i^2
P = [0,1,0,1,2,0,2,1,2,3,0,3,1,3,2,4,0,4,1,3,4,2,4,3,5,0,5,1,5,2,4,5,3,6,0,6,1,6,2,5,4,6,3,5];
Q = [0,0,1,1,0,2,1,2,2,0,3,1,3,2,3,0,4,1,4,3,2,4,3,4,0,5,1,5,2,5,4,3,5,0,6,1,6,2,6,4,5,3,6,5];

K = zeros(N,c);
for i = 1:N
    for j = 1:c
        K(i,j) = cos(pi*P(j)*Xn(i,1))*cos(pi*Q(j)*Xn(i,2));
    end
end

L = diag(P(1:c).^2 + Q(1:c).^2);
L = L.^1.0;
