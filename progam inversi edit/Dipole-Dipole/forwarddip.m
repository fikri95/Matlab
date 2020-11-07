function f=forwarddip(r,tebal,bvec)

suku1=forward(r,tebal,bvec);
turunanschl=turunans(r,tebal,bvec);
suku2=(bvec/2).*turunanschl;
hasil=suku1-suku2;


f=hasil;