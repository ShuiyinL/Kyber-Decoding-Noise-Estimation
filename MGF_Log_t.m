function [ MGF_log ] = MGF_Log_t(pmf_range,pmf_d,t)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
MGF_log=logdomain_sum_more(log(pmf_d)+t*pmf_range);
end

