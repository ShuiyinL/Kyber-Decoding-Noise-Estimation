function output= CBDistribution(k,m,n)

output=zeros(k,m);
for i=1:k
a=(randi(2,m,n)-1);
b=(randi(2,m,n)-1);
output(i,:)=sum(a-b,2)';
end

end

