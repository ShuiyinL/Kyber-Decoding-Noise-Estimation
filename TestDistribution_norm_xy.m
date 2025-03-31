clear all
close all

n=256;
weight=20;
eta=2;
num_samples=100000;


x_fixed=zeros(1,n);
x_fixed(1:weight)=1;


norm_sample_fixed=zeros(1,num_samples);
norm_sample_rand=zeros(1,num_samples);



for i=1:1:num_samples
    y=CBDistribution(1,n,eta); 
%   y=randn(1,n); %normal distribution
    
    x_rand=zeros(1,n);
    loc_ones=randperm(n,weight);
    x_rand(loc_ones)=1;
   
    
    norm_sample_fixed(i)=norm(x_fixed.*y);
    norm_sample_rand(i)=norm(x_rand.*y);
    
end


figure
h1=cdfplot(norm_sample_fixed);
set( h1, 'LineStyle', '-', 'Color', 'b','LineWidth',1);
hold on
h2=cdfplot(norm_sample_rand);
set( h2, 'LineStyle', '-', 'Color', 'r','LineWidth',1);
legend('CDF of |x_{fixed}*y|','CDF of |x_{rand}*y|','Location','best')
ylabel('CDF(z)');
xlabel('z');
title(['weight = ' num2str(weight)]);