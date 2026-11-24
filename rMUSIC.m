function [source_doa, Ra, Us] = rMUSIC(X,K)
% Optional outputs: actual sample covariance and its leading signal basis.

[M,N] = size(X);
R  = 1/N*X*X';
Ra = R;

[U,S] = svd(R);
if nargout > 2, Us = U(:,1:K); end
Un = U(:,K+1:end);
Gn = Un*Un';

a = zeros( 2*M-1, 1 )';

for i=-(M-1):(M-1) 
    a(i+M) = sum( diag(Gn,i) );
end

a1 = roots(a);

a2 = a1(abs(a1)<1);

[~, I] = sort( abs(abs(a2)-1) );

f = a2( I(1:K) );

for i = 1:K,
    source_doa(i) = -asin(angle(f(i))/pi) * 180/pi;
end

source_doa = sort(source_doa);
source_doa = source_doa(:);

end




