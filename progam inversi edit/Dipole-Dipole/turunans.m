function f=turunans(r,tebal,bvec)
r(1);
delx=0.001;
ourvec=zeros(size(bvec));
tesvec=zeros(size(bvec));
 tesvec1=zeros(size(bvec));
ourvec1=zeros(size(bvec));
dpds=zeros(size(bvec));
for i=1:length(bvec)
  B=bvec(i);
    B1=B+delx;
  ourvec(i)=testing(r,tebal,B);
  ourvec1(i)=testing(r,tebal,B1);
  tesvec(i)=r(1)+((B^2).*ourvec(i));
  tesvec1(i)= r(1)+((B1^2).*ourvec1(i));
    dpds(i)=tesvec1(i)-tesvec(i);
end

f=dpds;