function [mse,pr] = rmse(est, true)
% *********************************************************************
%input:
%DOA_hat: the estimated angles of signal
%DOA:the real angles of signal
%output:
%err:RMSE of the estimated errors
%sum_temp: the square of the err
%**********************************************************************
[K,Monte] = size(est);
sum_tmp=zeros(K,1);

est = sort(est,1);
true = sort(true,1);

for i=1:Monte,
    err_tmp = (est(:,i) - true').^2;
    sum_tmp = sum_tmp + err_tmp;
end
sum_tmp = sum_tmp/Monte;
mse = sum(sum_tmp);

t = 0;
true = true.'*ones(1,Monte);
err = abs(est - true);
for i = 1:Monte,
    tmp = find(err(:,i)<abs(true(1)-true(2))/2);
%     tmp = find(err(:,i)<1);
    if length(tmp) == K
        t = t + 1;
    end
end
pr = t/Monte;


end