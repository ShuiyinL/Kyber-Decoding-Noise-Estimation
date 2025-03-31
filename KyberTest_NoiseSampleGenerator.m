clear all
close all

%%%%%%%%%  Use to generate the Kyber decoding noise samples %%%%%%
% Reproduce the dataset used in the paper:
% Liu and A. Sakzad, “Lattice Codes for CRYSTALS-Kyber,” 2023.[Online]. Available: https://arxiv.org/abs/2308.13981
n=256;
q=3329;

flag=1; %1: Kyber512; 2: Kyber768, 3: Kyber1024

if flag==1
    k=2;
    eta1=3;
    eta2=2;
    d_tu=10;
    d_v=4;
%     d_tu=8;
elseif flag==2
    k=3;
    eta1=2;
    eta2=2;
    d_tu=10;
    d_v=4;
else
    k=4;
    eta1=2;
    eta2=2;
    d_tu=11;
    d_v=5;

end


%%%%%%%  Set d_v =12 to check Heuristic 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Heuristic 1: The major component of kyber decoding noise can be approximated by a Gaussian        
d_v =12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_itr=10000; %also determine the number of noise samples
cn_sample=zeros(1,num_itr);
noise_sample=zeros(1,num_itr);
noise_16norm=zeros(1,num_itr);
num_err=0;



for i=1:1:num_itr
    
%%%%%%%%%%%%%%%%%%%% Key generation: Public Key and Secrete Key %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matrix A
if k==2
    A_0={randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1};
elseif k==3
    A_0={randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1}; 
elseif k==4
    A_0={randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1;randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1,randi(q,1,n)-1};    
end

A_0T=A_0';
A=cell2mat(A_0);
A_T=cell2mat(A_0T);

%Secrete key s
s=CBDistribution(k,n,eta1);
s_T=reshape(s',1,[]);


% noise term e
e=CBDistribution(k,n,eta1);
e_T=reshape(e',1,[]);

%Public key t 
t_u = mod(PolyMatrixProduct(A,s,n,q)+e,q);
t_uT=reshape(t_u',1,[]);
t = t_u;

%%%%%%%%%%%%%%%%%%%%% Encryption %%%%%%%%%%%%%%%%%%%%%
%Plaintext
m=randi(2,1,n)-1;
% m=zeros(1,256); 

%r, e1, e2
t_T=reshape(t',1,[]);
r=CBDistribution(k,n,eta1);
e1=CBDistribution(k,n,eta2);
e2=CBDistribution(1,n,eta2);

%u: the first part of ciphertext
u_u = mod(PolyMatrixProduct(A_T,r,n,q)+e1,q);
u = CompressQ(u_u,q,d_tu);
% stu=PolyMatrixProduct(s_T,u_u,n,q);



%v: the second part ofciphertext
v_u=mod(PolyMatrixProduct(t_T,r,n,q)+e2+round(q/2)*m,q);
v = CompressQ(v_u,q,d_v);

%%%%%%%%%%%%%%%%%%%% Decryption  %%%%%%%%%%%%%%%%%%%%%
%Decompress (u,v)
u_d=DecompressQ(u,q,d_tu);
v_d=DecompressQ(v,q,d_v);
tempt=PolyMatrixProduct(s_T,u_d,n,q);

%Estimate the compression noise for u
cnoise=reshape(u_u-u_d,1,k*n);
cn_sample(i)=modFunc(cnoise(1,1),q);


%Estimate m
m_est=CompressQ(v_d-tempt,q,1);
num_err=num_err+sum(m_est~=m);

%%%%%%%%%%%%%%%%%%%%  Decoding noise samples %%%%%%%%%%%%%%%%%%%%%
decoding_noise=modFunc(v_d-tempt-round(q/2)*m,q);
noise_sample(i)=decoding_noise(1);
noise_16norm(i)=norm(decoding_noise(1:16));  %use for estimating the norms

end


%%%%%%%  Set d_v =12 at beginning to check Heuristic 1  %%%%%%
if d_v==12
    x_sample=(noise_sample-mean(noise_sample))/var(noise_sample)^.5;
    figure
    h1=cdfplot(x_sample);
    set( h1, 'LineStyle', '-', 'Color', 'b','LineWidth',1);
    hold on
    x_values = linspace(min(x_sample),max(x_sample));
    plot(x_values,normcdf(x_values,0,1),'r-','LineWidth',1)
    legend('Empirical CDF','Standard Normal CDF','Location','best')
    ylabel('CDF(x)');
    title('');
else
    var_noise=var(noise_sample/round(q/2)) % normalized decoding noise variance 
end

if d_tu==8 || d_tu==9
    var(cn_sample);
end

%%%% save the noise samples
% save kyberXXXX_16norm_n_XXXX.mat noise_16norm noise_sample




