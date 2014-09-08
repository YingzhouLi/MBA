function u = bfioHFFT(N,SL,EL,EXT,EPS,fun,f,mats,dir,stoplev)
% Hierarchical interpolation using fft2
if N<=2^stoplev
    % Direct evaluation
    % When N<64, bfio is not acurate
    kg = [-N/2:N/2-1];  %kg = [0:N-1];
    [k1,k2] = ndgrid(kg);
    ks = [k1(:)'; k2(:)'];
    xg = [0:1/N:(N-1)/N];
    [x1,x2] = ndgrid(xg);
    xs = [x1(:)'; x2(:)'];
    u = zeros(N^2,1);
    for cnt = 1:N^2
        u(cnt) = fun(N,xs(:,cnt),ks)*(f(:));
    end
    u = reshape(u,N,N);
else
    u = zeros(N,N);
    u = bfio_singlelevel_eval(N,SL,EL,EXT,EPS,fun,f,mats,dir);
    tmp = SpaceItp(N,N/2,bfioHFFT(N/2,SL-1,EL,EXT,EPS,fun,f(end/4+1:3*end/4,end/4+1:3*end/4),mats,dir,stoplev));
    if 0
        fprintf('compare direct method and bfio');
        checknum = 10;
        pos1 = floor(rand(1,checknum)*N)+1; pos2 = floor(rand(1,checknum)*N)+1;
        tmp(pos1,pos2)
        kg = [-N/4:N/4-1];
        [k1,k2] = ndgrid(kg);
        ks = [k1(:)'; k2(:)'];
        xg = [0:1/N:(N-1)/N];
        [x1,x2] = ndgrid(xg);
        xs = [x1(:)'; x2(:)'];
        ff = f(end/4+1:3*end/4,end/4+1:3*end/4);
        temp = zeros(N^2,1);
        for cnt = 1:N^2
            temp(cnt) = fun(N,xs(:,cnt),ks)*ff(:);
        end
        temp = reshape(temp,N,N);
        temp(pos1,pos2)
    end
    if 0
        fprintf('compare direct method and bfio');
        checknum = 30;
        app = zeros(1,checknum); 
        kg = [-N/4:N/4-1];
        [k1,k2] = ndgrid(kg);
        ks = [k1(:)'; k2(:)'];
        xg = [0:1/N:(N-1)/N];
        [x1,x2] = ndgrid(xg);
        xs = [x1(:)'; x2(:)'];
        ff = f(end/4+1:3*end/4,end/4+1:3*end/4);
        ext = zeros(1,checknum);
        for cnt = 1:checknum
            pos1 = floor(rand(1)*N); pos2 = floor(rand(1)*N);
            ext(cnt) = fun(N,[pos1/N;pos2/N],ks)*ff(:);
            app(cnt) = tmp(pos1+1,pos2+1); 
        end
        [app ;ext]
    end
    u = u + tmp;
end
