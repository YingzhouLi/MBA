function u = bfioDFFT(N,SL,EL,EPS,fun,f,mats,dir,stoplev)
% Direct interpolation using fft2
u = zeros(N,N);
NN = N;
lev = 1;
while N>2^stoplev
    if lev == 1
        u = bfio_eval(N,N,SL,EL,EPS,fun,f,mats,dir,dir);
    else
        tmp = SpaceItp(NN,N,bfio_eval(N,N,SL,EL,EPS,fun,f,mats,dir,dir));
        u = u + tmp;
    end
    N = N/2;
    f = f(end/4+1:3*end/4,end/4+1:3*end/4);
    SL = SL -1;
    lev = lev + 1;
end
% Direct evaluation
% When N<64, bfio is not acurate
kg = [-N/2:N/2-1];  %kg = [0:N-1];
[k1,k2] = ndgrid(kg);
ks = [k1(:)'; k2(:)'];
xg = [0:1/NN:(NN-1)/NN];
[x1,x2] = ndgrid(xg);
xs = [x1(:)'; x2(:)'];
tmp = zeros(NN^2,1);
for cnt = 1:NN^2
    tmp(cnt) = fun(NN,xs(:,cnt),ks)*(f(:));
end
u = u + reshape(tmp,NN,NN);





