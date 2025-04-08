clear all
close all

%The following code was used to address the review question:
%   *I invite the authors to experimentally observe the distribution norm $||x \cdot y||$ for random y when x has shape $(111 \ldots 000)$ vs $x$ is a random binary vector of the same weight
 
%We initially implemented the approach by considering the vectors (x, y) from \mathbb{Z}^n.  
%- There may have been a misunderstanding of the problem statement.  
%- The elements (x, y) might actually belong to the cyclotomic ring \mathbb{Z}[x] / (x^n + 1)."

%We would like to thank Prof. Leo Ducas for his comments and code on GitHub.  
%- We have incorporated the case of multiplication over the cyclotomic ring in Lines 49-50.  
%- The distribution of y has been selected according to Prof. Leo Ducas’s implementation.  
%- We are able to reproduce the Python plot from Prof. Leo Ducas's code on GitHub:
%  *see Norm_Plot_Python_Prof_Leo.pdf


n=256;
weight=20;
eta=2;
num_samples=1000;


x_fixed=zeros(1,n);
x_fixed(1:weight)=1;

norm_sample_fixed=zeros(1,num_samples);
norm_sample_rand=zeros(1,num_samples);


for i=1:1:num_samples
    y=x_fixed(randperm(n)); % code from Prof Leo Ducas  
%   y=CBDistribution(1,n,eta); Centred Binomial
%   y=randn(1,n); %normal distribution

    
    x_rand=zeros(1,n);
    loc_ones=randperm(n,weight);
    x_rand(loc_ones)=1;
   
    

    %Element-wise multiplication of (x, y) over \mathbb{Z}^n; review of the question (possible misunderstanding).
%     norm_sample_fixed(i)=norm(x_fixed.*y);
%     norm_sample_rand(i)=norm(x_rand.*y);

    %Multipilcation of (x,y) over the cyclotomic ring Z(x)/(x^n+1); 
    %Comments/Codes from Prof Leo Ducas
    norm_sample_fixed(i) = norm(cyclotomic_product(x_fixed,y,n));
    norm_sample_rand(i)  = norm(cyclotomic_product(x_rand,y,n));
    
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