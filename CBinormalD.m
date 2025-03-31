function [pmf_range,pmf_d_t] = CBinormalD( eta1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

pmf_d_t=zeros(1,2*eta1+1);
n=2*eta1;
for k=0:1:n
% nchoosek(2*eta1,1)
pmf_d_t(k+1)=nchoosek(n,k)*(1/2)^n;
pmf_range=-eta1:1:eta1;

end

