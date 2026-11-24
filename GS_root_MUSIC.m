function DOA = GS_root_MUSIC(X,P,true_DOA,resample_times,roll,pseudo_noise)
% Select the root-MUSIC estimate with the smallest DOA error.

    if nargin < 6
        pseudo_noise = Generate_pseudo_noise(X,P,resample_times,roll);
    end
    RMSE_min = Inf;

    for i = 0:resample_times
        if i == 0
            X_current = X;
        else
            X_current = X + pseudo_noise(:,:,i);
        end

        Estimated_DOA = rMUSIC(X_current,P);
        RMSE_current = rmse(Estimated_DOA,true_DOA).^0.5;
        if RMSE_current < RMSE_min
            RMSE_min = RMSE_current;
            DOA = Estimated_DOA;
        end
    end
end
