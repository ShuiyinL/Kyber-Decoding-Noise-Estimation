function [delta_BICM] = deltaBICM_com(n,t,delta, num_blocks)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t1=t+1;
log2e=log(2);
BICM_log=(delta+log2(num_blocks)-log2(n))*t1*log2e+gammaln(n+1)-gammaln(t1+1)-gammaln(n-t1+1);
%%https://stackoverflow.com/questions/24810597/how-to-calculate-a-large-combinatorial-function-in-matlab
for j=t1+1:1:min(t1+6,n)
    %tempt=(delta+4-log2(320))*j*log2e+gammaln(320+1)-gammaln(j+1)-gammaln(320-j+1);
    tempt=(delta+log2(num_blocks)-log2(n))*j*log2e+gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1);
    BICM_log=(logdomain_sum(BICM_log,tempt));
end

delta_BICM=log2(exp(1))*BICM_log;
end

