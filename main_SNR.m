function results = main_SNR(nT)
% Figure 2. Without arguments, run the 500 trials per point used in the paper.
if nargin == 0, nT = 500; end
oldPath = path;
oldRng = rng;
pathCleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
rngCleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
addpath(fileparts(mfilename('fullpath')), '-begin');
rng(1, 'twister');
M = 8;                              % Number of array elements
DOA = [5.5,16.2,30.7];
P = length(DOA);                    % Number of sources
K = 50;                            % Number of snapshots

SNR = linspace(-15,5,11);
beta = 0.5;
roll = 1;
resample_times = 50;

for iS = 1:length(SNR)

    snr = SNR(iS);

    for iT = 1:nT

        if rem(iT,nT/2) == 0
            fprintf('n = %d, Trials = %d, total = %d\n', ...
                iS, iT, (iS-1)*nT+iT);
        end

        X = StatSigGenerate(M, K, DOA, snr);
        pseudo_noise = Generate_pseudo_noise(X,P,resample_times,roll);
        [doa1(:,iT)] = GS_MUSIC(X,P,DOA,resample_times,roll,pseudo_noise);
        [doa2(:,iT)] = JWS_MUSIC(X,P,resample_times,beta,roll,pseudo_noise);
        [doa3(:,iT)] = GS_root_MUSIC(X,P,DOA,resample_times,roll,pseudo_noise);
        [doa4(:,iT)] = JWS_root_MUSIC(X,P,resample_times,beta,roll,pseudo_noise);
        [doa5(:,iT)] = MUSIC(X,P);
        [doa6(:,iT)] = rMUSIC(X,P);

    end

    RMSE1(iS) = rmse(doa1, DOA);
    RMSE2(iS) = rmse(doa2, DOA);
    RMSE3(iS) = rmse(doa3, DOA);
    RMSE4(iS) = rmse(doa4, DOA);
    RMSE5(iS) = rmse(doa5, DOA);
    RMSE6(iS) = rmse(doa6, DOA);

    % CRB
    [x1, A1, R_idl1, Rs1] = StatSigGenerate(M, K, DOA, snr);
    CRB(iS) = doa_crb_trace(A1,R_idl1,Rs1,DOA,K,1)*(180/pi)^2;

end

mz = 8;
lw = 2;

figure
semilogy(SNR, RMSE5.^0.5, '--v', 'markersize', mz, 'linewidth', lw); hold on;
semilogy(SNR, RMSE6.^0.5, '--p', 'markersize', mz, 'linewidth', lw);
semilogy(SNR, RMSE1.^0.5, '-d', 'markersize', mz, 'linewidth', lw);
semilogy(SNR, RMSE3.^0.5, '-^', 'markersize', mz, 'linewidth', lw);
semilogy(SNR, RMSE2.^0.5, '-o', 'markersize', mz, 'linewidth', lw);
semilogy(SNR, RMSE4.^0.5, '-*', 'markersize', mz, 'linewidth', lw);
semilogy(SNR, CRB.^0.5, 'k-', 'linewidth', lw);
xlabel('SNR(dB)'); ylabel('RMSE (degree)');
legend('MUSIC','Root-MUSIC','GS-MUSIC','GS-root-MUSIC', ...
    'JWS-MUSIC','JWS-root-MUSIC','CRB','Location','best');
grid on;
results = table(SNR(:),sqrt(RMSE5(:)),sqrt(RMSE6(:)),sqrt(RMSE1(:)), ...
    sqrt(RMSE3(:)),sqrt(RMSE2(:)),sqrt(RMSE4(:)),sqrt(CRB(:)), ...
    'VariableNames',{'SNR_dB','MUSIC','RootMUSIC','GS_MUSIC','GS_RootMUSIC', ...
    'JWS_MUSIC','JWS_RootMUSIC','CRB'});
end

