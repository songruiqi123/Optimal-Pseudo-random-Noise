function varianceSum = doa_crb_trace(A, R, sourceCov, anglesDeg, snapshots, noiseVariance)
% Trace of the stochastic DOA CRB for the half-wavelength ULA, in rad^2.
% Stoica and Nehorai (1989), doi:10.1109/29.17564.
% A scalar sourceCov denotes equal source powers (sourceCov times identity).
sensor = (0:size(A,1)-1).';
derivative = 1i*pi*sensor .* A .* cosd(anglesDeg(:).');
perpendicular = derivative - A*(A\derivative);
signal = A*sourceCov;
signalMetric = signal'*(R\signal);
information = (2*snapshots/noiseVariance) * ...
    real((perpendicular'*perpendicular) .* conj(signalMetric));
factor = chol((information+information.')/2);
varianceSum = norm(factor\eye(numel(anglesDeg)),'fro')^2;
end
