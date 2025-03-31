function C = PolyMatrixProduct(B,e,n,q)
[r,c]=size(B);
k=c/n;
C=zeros(r,n);
Xn_1=[1, zeros(1,n-1), 1];
for i=1:1:r
    tempt=0;
    for j=1:1:k
        [~,r_t] = deconv(conv(B(i,(j-1)*n+1:j*n),e(j,:)),Xn_1);
        tempt=tempt+mod(r_t(end-n+1:1:end),q);
    end
    
    C(i,:)=mod(tempt,q);
end
end

