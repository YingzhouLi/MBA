function u = bfioChebyshev(N,SL,EL,EXT,EPS,fun,f,mats,dir,dirlev,stoplev,lev)
%N
NN = N*2^(lev-1);
if N<=2^stoplev
    % Direct evaluation
    % When N<64, bfio is not acurate
    kg = [-N/2:N/2-1];  %kg = [0:N-1];
    [k1,k2] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'];
    xg = [0:1/NN:(NN-1)/NN];
    [x1,x2] = ndgrid(xg);
    xs = [x1(:)'; x2(:)'];
    u = zeros(NN^2,1);
    for cnt = 1:NN^2
        u(cnt) = fun(NN,xs(:,cnt),ks)*(f(:));
    end
    u = reshape(u,NN,NN);
else
    u = zeros(NN,NN);
    u = bfio_multilevel_eval(N,SL,EL,EXT,EPS,fun,f,mats,dir,dirlev,lev);
    %u = u + SpaceItp(N/2,itp,bfioNew(N/2,SL-1,EL,EXT,EPS,fun,f(end/4+1:3*end/4,end/4+1:3*end/4),mats,dir));
    tmp = bfioChebyshev(N/2,SL-1,EL,EXT,EPS,fun,f(end/4+1:3*end/4,end/4+1:3*end/4),mats,dir,dirlev,stoplev,lev+1);
    if 0
        fprintf('compare direct method and bfio');
        lev
        checknum = 10;
        pos1 = floor(rand(1,checknum)*NN)+1; pos2 = floor(rand(1,checknum)*NN)+1;
        tmp(pos1,pos2)
        kg = [-N/4:N/4-1];
        [k1,k2] = ndgrid(kg);
        ks = [k1(:)'; k2(:)'];
        xg = [0:1/NN:(NN-1)/NN];
        [x1,x2] = ndgrid(xg);
        xs = [x1(:)'; x2(:)'];
        ff = f(end/4+1:3*end/4,end/4+1:3*end/4);
        temp = zeros(NN^2,1);
        for cnt = 1:NN^2
            temp(cnt) = fun(N,xs(:,cnt),ks)*ff(:);
        end
        temp = reshape(temp,NN,NN);
        temp(pos1,pos2)
    end
    if 1
        fprintf('compare direct method and bfio');
        lev
        checknum = 30;
        app = zeros(1,checknum); 
        kg = [-N/4:N/4-1];
        [k1,k2] = ndgrid(kg);
        ks = [k1(:)'; k2(:)'];
        xg = [0:1/NN:(NN-1)/NN];
        [x1,x2] = ndgrid(xg);
        xs = [x1(:)'; x2(:)'];
        ff = f(end/4+1:3*end/4,end/4+1:3*end/4);
        ext = zeros(1,checknum);
        for cnt = 1:checknum
            pos1 = floor(rand(1)*NN); pos2 = floor(rand(1)*NN);
            ext(cnt) = fun(N,[pos1/NN;pos2/NN],ks)*ff(:);
            app(cnt) = tmp(pos1+1,pos2+1); 
        end
        [app ;ext]
    end
    u = u + tmp;
end
