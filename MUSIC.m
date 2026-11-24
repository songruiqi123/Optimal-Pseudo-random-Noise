function [source_doa, Ra, Us] = MUSIC(X,K)
% Optional outputs: actual sample covariance and its leading signal basis.

[M,N] = size(X);
R  = 1/N*X*X';
Ra = R;

[U,S] = svd(R);
if nargout > 2, Us = U(:,1:K); end
Un = U(:,K+1:end);
% MUSIC spatial spectrum
search_angle = -90:0.1:90;

A = exp(1j*pi*(0:M-1)'*sind(search_angle));

P_MUSIC = 1./max(sum(abs(Un'*A).^2,1),eps);

% Peak detection
[~,locs] = findpeaks(P_MUSIC);

[~,I] = sort(P_MUSIC(locs),'descend');

% Select up to K detected peaks to avoid indexing past the available peaks.
num_selected = min(K, numel(locs));
source_doa = search_angle(locs(I(1:num_selected)));

source_doa = sort(source_doa);
% If fewer than K peaks are detected, pad the DOA estimates with 90 degrees.
source_doa(end+1:K) = 90;
source_doa = source_doa(:);
end




