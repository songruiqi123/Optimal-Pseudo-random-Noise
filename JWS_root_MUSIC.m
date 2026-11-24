function DOA = JWS_root_MUSIC(X,P,resample_times,beta,roll,pseudo_noise)
% Select the DOA estimate with the lowest joint weighted cost.

    [M,N] = size(X);
    R = X*X'/N;
    [~,V] = svd(R);
    eigenvalues = real(diag(V));
    noise_power = mean(eigenvalues(P+1:end));

    Estimated_DOA = nan(P,resample_times+1);
    cost_function = inf(1,resample_times+1);

    for i = 1:resample_times+1
        if i == 1
            X_current = X;
        elseif nargin >= 6
            X_current = X + pseudo_noise(:,:,i-1);
        else
            Y_k = sqrt(roll*noise_power) * (randn(M,N)+1j*randn(M,N))/sqrt(2);
            X_current = X + Y_k;
        end
        DOA_current = rMUSIC(X_current,P);
        Estimated_DOA(:,i) = DOA_current;

        A = exp(1j*pi*(0:M-1)'*sind(DOA_current(:).'));
        [Q,~] = qr(A,0);
        P_signal = Q*Q';
        P_noise = eye(M)-P_signal;

        % Cross-correlation cost from the original covariance matrix.
        T = P_signal*R*P_noise;
        power_signal = real(trace(R*P_signal))-P*noise_power;
        cross_cost = norm((T+T')/power_signal,'fro')^2;

        % Weighted subspace fitting cost from the candidate covariance matrix.
        RR = X_current*X_current'/N;
        [U,V] = svd(RR);
        eigenvalues = real(diag(V));
        signal_eigenvalues = eigenvalues(1:P);
        gaps = signal_eigenvalues-mean(eigenvalues(P+1:end));
        weights = gaps.^2./signal_eigenvalues;
        U_signal = U(:,1:P);
        W = diag(weights/sum(weights));
        wsf_cost = real(trace(P_noise*U_signal*W*U_signal'));
        cost_function(i) = beta*cross_cost+(1-beta)*wsf_cost;
    end

    valid = isfinite(cost_function);
    cost_function(~valid) = Inf;
    Estimated_DOA(:,~valid) = NaN;
    [~,idx] = min(cost_function);
    DOA = Estimated_DOA(:,idx);
end
