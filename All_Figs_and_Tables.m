clear all
close all

%%%%%%%% Dataset and codes used to generate All figures and tables in the paper:
% Liu and A. Sakzad, “Lattice Codes for CRYSTALS-Kyber,” 2023.
% Table 8 summarizes existing results from the literature, thus is not included.
% The calculations of CER are straightforward, and thus omitted.


q=3329;
n=256;
scalar=round(q/2);

flag=1; %1: Kyber512; 2: Kyber768, 3: Kyber1024

if flag==1     %1: Kyber512
    
    fprintf('Summary of DFRs of Kyber512 using different encoders \n')
    
    eta1=3;
    eta2=2;
    k=2;
    
    d=4;
    b_q=round(q/2^(d+1));
    
    pmf_compress =[0.0385    0.3082    0.3095    0.3062    0.0376];
    pmf_compress_range =[-2    -1     0     1     2];
    
    var_cu=pmf_compress_range.^2*pmf_compress'-(pmf_compress_range*pmf_compress')^2;
    sigma_s=(((eta2/2+round(var_cu,1))*eta1/2+eta1*eta1/4)*k*n+eta2/2)^.5;  %round var_cu according to [7], also see Table 6
    
    %%%%%%%%%%%% Compute F(t) using the Python code from [25] https://github.com/pq-crystals/security-estimates
    %Step 1: one change in the codes:  set d_v=12; 
    %        *in Kyber.py:  replace KyberParameterSet(256, 2, 3, 3, 3329, 2**12, 2**10, 2**4, ke_ct=2)
    %                       to      KyberParameterSet(256, 2, 3, 3, 3329, 2**12, 2**10, 2**12, ke_ct=2)
    %        *run Kyber.py, obtain the DFR
    %Step 2: F(t) = DFR/2^8
    %        *DFR was computed using the union bound: DFR=F(t)*2^8
    %        *F(t): DFR for a single information bit
    
    ep_sim=-180; %F(t)=2^(-180)
    
    theta=0.2490; %optimal theta in Remark 2
    
    
elseif flag==2  %2: Kyber768
    
    fprintf('Summary of DFRs of Kyber768 using different encoders \n')
    
    eta1=2;
    eta2=2;
    k=3;
    
    
    d=4;
    b_q=round(q/2^(d+1));
    pmf_compress =[0.0385    0.3082    0.3095    0.3062    0.0376];
    pmf_compress_range =[-2    -1     0     1     2];
    
    var_cu=pmf_compress_range.^2*pmf_compress'-(pmf_compress_range*pmf_compress')^2;
    
    sigma_s=(((eta2/2+round(var_cu,1))*eta1/2+eta1*eta1/4)*k*n+eta2/2)^.5;  %round var_cu according to [7], also see Table 6
    
    %%%%%%%%%%%% Compute F(t) using the Python code from [25] https://github.com/pq-crystals/security-estimates
    %Step 1: one change in the codes:  set d_v=12; 
    %        *in Kyber.py:  replace KyberParameterSet(256, 3, 2, 2, 3329, 2**12, 2**10, 2**4)
    %                       to      KyberParameterSet(256, 3, 2, 2, 3329, 2**12, 2**10, 2**12)
    %        *run Kyber.py, obtain the DFR
    %Step 2: F(t) = DFR/2^8
    %        *DFR was computed using the union bound: DFR=F(t)*2^8
    %        *F(t): DFR for a single information bit
    
    ep_sim=-214; %F(t)=2^(-214)
    
    theta=0.2990; %optimal theta in Remark 2
    
    
elseif flag==3   %3: Kyber1024
    
    fprintf('Summary of DFRs of Kyber1024 using different encoders \n')

    eta1=2;
    eta2=2;
    k=4;
    
    d=5;
    b_q=round(q/2^(d+1));
    pmf_compress =[0.1916    0.6161    0.1923];
    pmf_compress_range =[ -1     0     1];
    
    var_cu=pmf_compress_range.^2*pmf_compress'-(pmf_compress_range*pmf_compress')^2;    
    sigma_s=(((eta2/2+round(var_cu,2))*eta1/2+eta1*eta1/4)*k*n+eta2/2)^.5;   %round var_cu, also see Table 6
    
    %%%%%%%%%%%% Compute F(t) using the Python code from [25] https://github.com/pq-crystals/security-estimates
    %Step 1: one change in the codes:  set d_v=12; 
    %        *in Kyber.py:  replace KyberParameterSet(256, 4, 2, 2, 3329, 2**12, 2**11, 2**5)
    %                       to      KyberParameterSet(256, 4, 2, 2, 3329, 2**12, 2**11, 2**12)
    %        *run Kyber.py, obtain the DFR
    %Step 2: F(t) = DFR/2^8
    %        *DFR was computed using the union bound: DFR=F(t)*2^8
    %        *F(t): DFR for a single information bit
    
    ep_sim=-201; %F(t)=2^(-201)
    
    theta=0.3; %optimal theta in Remark 2
    
    
    
    %%%%%% Figure 2 and Figure 3   
    load kyber1024_16norm_n_1000000.mat noise_16norm noise_sample  
    range_dist=floor(min(noise_16norm)):10: ceil(max(noise_16norm));
    j=1;
    l=16;
    
    sim_l_bound=zeros(1,length(range_dist));
    CLT_l_bound=zeros(1,length(range_dist));
    for i=range_dist
        CLT_l_bound(j)=(marcumq(l^.5*b_q/sigma_s,i/sigma_s,ceil(l/2)));
        sim_l_bound(j)=(sum(noise_16norm>i)/length(noise_16norm));
        j=1+j;
    end
    
    %Figure 2
    figure
    x_sample=(noise_sample-mean(noise_sample))/var(noise_sample)^.5;
    h1=cdfplot(x_sample);
    set( h1, 'LineStyle', '-', 'Color', 'b','LineWidth',1);
    hold on
    x_values = linspace(min(x_sample),max(x_sample));
    plot(x_values,normcdf(x_values,0,1),'r-','LineWidth',1)
    legend('Empirical CDF','Standard Normal CDF','Location','best')
    ylabel('CDF(x)');
    title('');
    
    %Figure 3
    figure
    plot(range_dist,log2(CLT_l_bound),'-b',range_dist,log2(sim_l_bound),'-r','linewidth',1.5);
    grid on
    xlabel('$x$','Interpreter','latex')
    ylabel('$\log_2(\Pr(\|n_e^{(16)}\|>x))$','Interpreter','latex')
    legend('Upper Bound in Lemma 4','Empirical Tail Distribution')
    
end


%%%%%% Figure 1 and Table 6
load kyber512_n_100000.mat
kyber512_Normalized_Noise_Variance=round(var(noise_sample/scalar),4);

load kyber768_n_100000.mat 
kyber768_Normalized_Noise_Variance=round(var(noise_sample/scalar),4);
%Figure 1
figure
normplot(noise_sample)

load kyber1024_n_100000.mat 
kyber1024_Normalized_Noise_Variance=round(var(noise_sample/scalar),4);



%%%%%% Table 2 when flag > 1
[pmf_range_1,pmf_d_1] = CBinormalD(eta1);
[pmf_range_2,pmf_d_2] = CBinormalD(eta2);
%Distribution of T1 for KYBER768 and KYBER1024
[pmf_range_T1,pmf_d_T1] =ProdPMF( pmf_range_1,pmf_d_1,pmf_range_1,pmf_d_1 );


%%%%%% Table 3 when flag ==3
[pmf_range_sum,pmf_d_sum] = SumPMF(pmf_range_2,pmf_d_2, pmf_compress_range,pmf_compress);
%Distribution of T2 for KYBER1024
[pmf_range_T2,pmf_d_T2] =ProdPMF( pmf_range_1,pmf_d_1,pmf_range_sum,pmf_d_sum );


%%%%%% Table 4
%Worst-case DFR
kn=k*n;
pmf_range_T3=[-b_q:1:b_q];
pmf_d_T3=1/(2*b_q+1)*ones(1,2*b_q+1);

%Lemma 2 (Worst-case upper bound on DFR).
M_T1=MGF_Log_t(pmf_range_T1,pmf_d_T1,theta)*kn;
M_T2=MGF_Log_t(pmf_range_T2,pmf_d_T2,theta)*kn;
M_T4=MGF_Log_t(pmf_range_2,pmf_d_2,theta);
M_T3=MGF_Log_t(pmf_range_T3,pmf_d_T3,theta);
DFR_WorstCase= log2(2*n)+log2(exp(1))*(M_T1+M_T2+M_T3+M_T4-theta*round(q/4))

%Corollary 1 (Average-case upper bound on DFR)
DFR_AverageCase=log2(2*n*qfunc((round(q/4)-round(q/2^(d+1)))/sigma_s)) 

%%%%%% Table 5
est_normal_tail=log2(qfunc(round(q/4)/(k*n*eta1^2/4+k*n*eta1/2*(eta2/2+var_cu)+eta2/2)^.5));
Kolmogorov_distance=logdomain_diff(ep_sim,est_normal_tail) 


%%%%%% Table 7
%BW16
ell=16;
decoding_radius=0.7067*round(q/2);%[31]
DFR_BW16=log2(n/ell)+log2(marcumq(ell^.5*round(q/2^(d+1))/sigma_s,decoding_radius/sigma_s,ceil(ell/2)))

%Leech
ell=24;
decoding_radius=0.7067*round(q/2);%[31]
% DFR of 11 Leech blocks > DFR of 10 Leech + 1 BW16
DFR_Leech24=log2(ceil(n/ell))+log2(marcumq(ell^.5*round(q/2^(d+1))/sigma_s,decoding_radius/sigma_s,ceil(ell/2)))  


%%%%%% Table 9: choose the value of d_u: 9 or 8
d_u=8;

if d_u==9
    var_cu=3.8; 
else
    var_cu=14.1;
end
sigma_s=(((eta2/2+var_cu)*eta1/2+eta1*eta1/4)*k*n+eta2/2)^.5;  
ell=16;
decoding_radius=0.7067*round(q/2); %[31]
DFR_BW16_one_block=log2(marcumq(ell^.5*round(q/2^(d+1))/sigma_s,decoding_radius/sigma_s,ceil(ell/2)));

if d_u==9
    DFR_BICM_d_u_9 = deltaBICM_com(320,7,DFR_BW16_one_block, n/ell)
else
    DFR_BICM_d_u_8 = deltaBICM_com(320,7,DFR_BW16_one_block, n/ell)
end



