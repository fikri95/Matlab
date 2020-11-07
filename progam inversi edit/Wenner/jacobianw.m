%Function Jacobian%
function[A]=jacobian(r,tebal,lr,lt,ds,r1,bvec)
par=0.1;
r2=r;
for i2=1:lr
    r2(i2)=r(i2)*(par)+r(i2);
    ro2=forwardw(r2,tebal,bvec);
    A1(:,i2)=[(ro2-r1)/(r(i2)*(par))]*r(i2)./ds;
    r2=r;
end
t2=tebal;
for i3=1:lt
    t2(i3)=tebal(i3)*(par)+tebal(i3);
    ro3=forwardw(r,t2,bvec);
    A2(:,i3)=[(ro3-r1)/(tebal(i3)*(par))]*tebal(i3)./ds;
    t2=tebal;
end
A=[A1 A2];
return