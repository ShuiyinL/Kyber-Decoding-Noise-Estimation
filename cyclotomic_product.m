function C = cyclotomic_product(a,b,n)
%Comments/Codes from Prof Leo Ducas
c=zeros(1,2*n);

for i=1:1:n  
    c(i:i+n-1)=a(i)*b+c(i:i+n-1);
end

C=c(1:n)-c(n+1:end);

end

