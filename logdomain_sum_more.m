function T=logdomain_sum_more(x)

it=max(size(x));


if isempty(x)==1
     T=0;
elseif it==1
    T=x;
    
else
T=logdomain_sum(x(1),x(2));
for i=3:1:it
    
    T=logdomain_sum(x(i),T);
    
end
end