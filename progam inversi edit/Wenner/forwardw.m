function f=forwardw(r,t,bvec)
a=(2/3).*bvec;
a2=2*a;
for i=1:length(bvec)
  B=a(i);
  BB=a2(i);
  ourvec1(i)=testingw1(r,t,B);
  ourvec2(i)=testingw2(r,t,BB);
  truevec(i)=BB.*(ourvec1(i)-ourvec2(i));
end
f=truevec;