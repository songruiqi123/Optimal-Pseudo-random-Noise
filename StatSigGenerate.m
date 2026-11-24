function [x, A, R_idl, Rs] = StatSigGenerate(M, N, DOA, SNR)

jj = sqrt(-1);

K = length(DOA);

sigma_s2 = 10.^(SNR/10);

sigma_n2 = 1;

nt = sqrt(0.5) * (randn(M,N) + jj*randn(M,N));

nt = nt/sqrt(trace(nt*nt'/N)/M);

st = sqrt(0.5) * (randn(K,N) + jj*randn(K,N));
st = st/sqrt(trace(st*st'/N)/K);

st = diag(sqrt(sigma_s2)) * st;

A = exp(jj*pi*[0:M-1]'*sind(DOA));

x = A*st + nt*sqrt(sigma_n2);

Rs = diag(sigma_s2);

R_idl = A* Rs * A' + sigma_n2 * eye(M);



