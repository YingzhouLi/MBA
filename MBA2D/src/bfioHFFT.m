function u = bfioHFFT(Nx,Nk,SL,EL,EPS,fun,f,mats,dir,stoplev)
% Hierarchical interpolation using fft2
if Nk<=2^stoplev
    % Direct evaluation
    % When N<64, bfio is not acurate
    kg = [-Nk/2:Nk/2-1];  %kg = [0:N-1];
    [k1,k2] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'];
    xg = [0:1/Nk:(Nk-1)/Nk];
    [x1,x2] = ndgrid(xg);
    xs = [x1(:)'; x2(:)'];
    u = zeros(Nk^2,1);
    for cnt = 1:Nk^2
        u(cnt) = fun(xs(:,cnt),ks)*(f(:));
    end
    u = reshape(u,Nk,Nk);
else
    u = bfio_eval(Nk,Nk,SL,EL,EPS,fun,f,mats,dir,dir);
    tmp = SpaceItp(Nk,Nk/2,bfioHFFT(Nx,Nk/2,SL-1,EL,EPS,fun,f(end/4+1:3*end/4,end/4+1:3*end/4),mats,dir,stoplev));
    u = u + tmp;
end
