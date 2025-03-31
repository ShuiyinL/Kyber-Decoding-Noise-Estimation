function [pmf_range,pmf_d] = SumPMF(pmf_range_1,pmf_d_1, pmf_range_2,pmf_d_2 )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

eta1=max(pmf_range_1);
eta2=max(pmf_range_2);
bd=eta1+eta2;
bd_left=min(pmf_range_1)+min(pmf_range_2);

% pmf_d=zeros(1,2*bd+1);
% pmf_range=-bd:1:bd;
pmf_d=zeros(1,bd-bd_left+1);
pmf_range=bd_left:1:bd;

sum_d_tempt=[];
sum_range_tempt=[];

% for i=1:1:2*eta2+1
for i=1:1:length(pmf_range_2)
    sum_range_tempt=[sum_range_tempt,pmf_range_1+pmf_range_2(i)];
    sum_d_tempt=[sum_d_tempt,pmf_d_1*pmf_d_2(i)];
end

% for ii=1:1:2*bd+1
for ii=1:1:bd-bd_left+1
%     pmf_d(ii)=sum(x==pmf_range(ii))/num_sample;
   pmf_d(ii)= sum(sum_d_tempt((sum_range_tempt==pmf_range(ii))));
    
end
pmf_range(pmf_d==0)=[];
pmf_d(pmf_d==0)=[];
pmf_d=pmf_d/sum(pmf_d);


end

