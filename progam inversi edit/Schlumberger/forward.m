function f=forward(r,tebal,bvec)
r(1);
for i=1:length(bvec),
  B=bvec(i);
  ourvec(i)=testing(r,tebal,B);
  tesvec1(i)= r(1)+((B^2).*ourvec(i));
end

f=tesvec1;