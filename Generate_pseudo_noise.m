function pseudo_noise = Generate_pseudo_noise(X,P,resample_times,roll)
% Generate the pseudo-random noise shared by all methods in one Monte Carlo trial.

    [M,N] = size(X);
    R = X*X'/N;
    power = svd(R);
    noise_power = mean(power(P+1:end));

    pseudo_noise = zeros(M,N,resample_times);
    for i = 1:resample_times
        pseudo_noise(:,:,i) = sqrt(roll*noise_power) * ...
            (randn(M,N) + 1j*randn(M,N))/sqrt(2);
    end

end
