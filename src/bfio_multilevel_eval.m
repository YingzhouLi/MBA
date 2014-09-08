function u = bfio_eval(N,SL,EL,EXT,EPS,fun,f,mats,dir,dirlev,lev)
if EL<2
    error('EL has to be >= 2');
end
NN = N*2^(lev-1);
%NN in x and N in k
NG = EPS;

grid = bfio_grid(EPS); %0-1 grid

u = zeros(NN,NN);

ML = floor((SL+EL)/2);

nz = 2^EL; %number of zones
zB = N/nz;

for z1=0:nz-1
    for z2=0:nz-1
        if z1<nz/4 | z2 < nz/4 | z1 > 3*nz/4-1 | z2 > 3*nz/4-1
            %fprintf(1,'%d %d\n', z1,z2);
            k1stt = z1*zB;      k1end = (z1+1)*zB;
            k2stt = z2*zB;      k2end = (z2+1)*zB;
            
            NOW = [];
            for ell=SL:-1:ML
                nk = 2^ell;    nx = N/nk;%nk: number of boxes in k; nx: number of boxes in x; nk*nx = N;
                kB = N/nk;    xB = 1/nx; %kB: size of boxes in k; xB: size of boxes in x; kB*xB = 1;
                %LLLLLLL
                PRE = NOW;
                NOW = cell(nk,nk);
                for k1=(k1stt/kB):(k1end/kB-1)
                    for k2=(k2stt/kB):(k2end/kB-1)
                        NOW{k1+1,k2+1} = cell(nx,nx);
                        for x1=0:nx-1
                            for x2=0:nx-1
                                NOW{k1+1,k2+1}{x1+1,x2+1} = zeros(NG,NG);
                            end
                        end
                    end
                end
                %----
                %in unit box
                [so1,so2] = ndgrid([0:kB-1]/kB);
                [ko1,ko2] = ndgrid(grid);        %[co1,co2] = ndgrid(grid*kB/2-N/2);
                
                for k1=(k1stt/kB):(k1end/kB-1)
                    for k2=(k2stt/kB):(k2end/kB-1)
                        kc1 = (k1+1/2)*kB - N/2;
                        kc2 = (k2+1/2)*kB - N/2;
                        
                        for x1=0:nx-1
                            for x2=0:nx-1
                                xc1 = (x1+1/2)*xB;
                                xc2 = (x2+1/2)*xB;
                                %--------------
                                if(ell==SL)
                                    s1 = [0:kB-1] + k1*kB - N/2;
                                    s2 = [0:kB-1] + k2*kB - N/2;
                                    %get
                                    all = f(s1+N/2+1,s2+N/2+1);
                                    %scale
                                    trg = [xc1;xc2];
                                    sp1 = (so1+k1)*kB - N/2;
                                    sp2 = (so2+k2)*kB - N/2;
                                    src = [sp1(:)'; sp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,kB,kB);
                                    all = all.*scl;
                                    %transform
                                    all = dir*all*dir';
                                    %scale
                                    trg = [xc1;xc2];
                                    kp1 = (ko1+k1)*kB - N/2;
                                    kp2 = (ko2+k2)*kB - N/2;
                                    src = [kp1(:)'; kp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,NG,NG);
                                    all = all./scl;
                                    %put
                                    %fprintf(1,'%d %d %d %d %d %d\n', k1,k2,x1,x2,real(all(1,1)),imag(all(1,1)));
                                    NOW{k1+1,k2+1}{x1+1,x2+1} = all;
                                    if(0)
                                        src = [sp1(:)'; sp2(:)'];
                                        tmp0 = fun(N,trg,src);                tmp0 = reshape(tmp0,kB,kB);
                                        tmp1 = fun(N,trg-3/N,src);                tmp1 = reshape(tmp1,kB,kB);
                                        fprintf(1,'%d %d %d %d\n', k1,k2,x1,x2);
                                        imagesc(real(tmp1./tmp0)); colorbar; pause(0.2);
                                    end
                                else
                                    p1 = floor(x1/2);              p2 = floor(x2/2);
                                    all = zeros(NG,NG);
                                    for a=0:1
                                        for b=0:1
                                            c1 = 2*k1+a;                  c2 = 2*k2+b;
                                            %get
                                            tmp = PRE{c1+1,c2+1}{p1+1,p2+1};
                                            %scale
                                            trg = [xc1;xc2];
                                            cp1 = (ko1+c1)*kB/2 - N/2;
                                            cp2 = (ko2+c2)*kB/2 - N/2;
                                            src = [cp1(:)'; cp2(:)'];
                                            scl = fun(N,trg,src);                  scl = reshape(scl,NG,NG);
                                            tmp = tmp.*scl;
                                            %transform
                                            all = all + mats{a+1}*tmp*mats{b+1}';
                                        end
                                    end
                                    %scale
                                    trg = [xc1;xc2];
                                    kp1 = (ko1+k1)*kB - N/2;
                                    kp2 = (ko2+k2)*kB - N/2;
                                    src = [kp1(:)'; kp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,NG,NG);
                                    all = all./scl;
                                    %put
                                    %fprintf(1,'%d %d %d %d %d %d\n', k1,k2,x1,x2,real(all(1,1)),imag(all(1,1)));error('w');
                                    NOW{k1+1,k2+1}{x1+1,x2+1} = all;
                                end
                            end
                        end
                        %x
                    end
                end
                %k
                clear PRE;
            end
            clear c1 c2 p1 p2;
            %--------------------------------
            if(1)
                ell = ML;
                nk = 2^ell;    nx = N/nk;
                kB = N/nk;    xB = 1/nx;
                
                %in unit box
                [ko1,ko2] = ndgrid(grid);
                [xo1,xo2] = ndgrid(grid);
                
                for k1=(k1stt/kB):(k1end/kB-1)
                    for k2=(k2stt/kB):(k2end/kB-1)
                        kc1 = (k1+1/2)*kB - N/2;
                        kc2 = (k2+1/2)*kB - N/2;
                        for x1=0:nx-1
                            for x2=0:nx-1
                                xc1 = (x1+1/2)*xB;
                                xc2 = (x2+1/2)*xB;
                                %switch
                                kp1 = (ko1+k1)*kB - N/2;
                                kp2 = (ko2+k2)*kB - N/2;
                                src = [kp1(:)'; kp2(:)'];
                                xp1 = (xo1+x1)*xB;
                                xp2 = (xo2+x2)*xB;
                                trg = [xp1(:)'; xp2(:)'];
                                all = fun(N,trg,src) * NOW{k1+1,k2+1}{x1+1,x2+1}(:);
                                NOW{k1+1,k2+1}{x1+1,x2+1} = reshape(all, NG,NG);
                                %fprintf(1,'%d %d %d %d %d %d\n', k1,k2,x1,x2,real(all(1,1)),imag(all(1,1)));
                            end
                        end
                    end
                end
            end
            %--------------------------------
            for ell=ML:-1:EL
                nk = 2^ell;    nx = N/nk;
                kB = N/nk;    xB = 1/nx;
                %LLLLLLLLL
                NXT = cell(nk/2,nk/2);
                for k1=(k1stt/(2*kB)):(k1end/(2*kB)-1)
                    for k2=(k2stt/(2*kB)):(k2end/(2*kB)-1)
                        NXT{k1+1,k2+1} = cell(2*nx,2*nx);
                        for x1=0:2*nx-1
                            for x2=0:2*nx-1
                                NXT{k1+1,k2+1}{x1+1,x2+1} = zeros(NG,NG);
                            end
                        end
                    end
                end
                %--
                
                NT = 2^(EL+lev-1);
                [to1,to2] = ndgrid([0:NT-1]/NT);
                [xo1,xo2] = ndgrid(grid);
                
                for k1=(k1stt/kB):(k1end/kB-1)
                    for k2=(k2stt/kB):(k2end/kB-1)
                        kc1 = (k1+1/2)*kB - N/2;
                        kc2 = (k2+1/2)*kB - N/2;
                        for x1=0:nx-1
                            for x2=0:nx-1
                                xc1 = (x1+1/2)*xB;
                                xc2 = (x2+1/2)*xB;
                                %--------------
                                if(ell~=EL)
                                    all = NOW{k1+1,k2+1}{x1+1,x2+1};
                                    %scale
                                    src = [kc1; kc2];
                                    xp1 = (xo1+x1)*xB;
                                    xp2 = (xo2+x2)*xB;
                                    trg = [xp1(:)'; xp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,NG,NG);
                                    all = all./scl;
                                    %--
                                    q1 = floor(k1/2);            q2 = floor(k2/2);
                                    for a=0:1
                                        for b=0:1
                                            d1 = 2*x1+a;                  d2 = 2*x2+b;
                                            %transform
                                            tmp = mats{a+1}'*all*mats{b+1};
                                            %scale
                                            src = [kc1; kc2];
                                            dp1 = (xo1+d1)*xB/2;
                                            dp2 = (xo2+d2)*xB/2;
                                            trg = [dp1(:)'; dp2(:)'];
                                            scl = fun(N,trg,src);              scl = reshape(scl,NG,NG);
                                            tmp = tmp.*scl;
                                            %put
                                            NXT{q1+1,q2+1}{d1+1,d2+1} = NXT{q1+1,q2+1}{d1+1,d2+1} + tmp;
                                        end
                                    end
                                    %ab
                                else
                                    all = NOW{k1+1,k2+1}{x1+1,x2+1};
                                    %scale
                                    src = [kc1; kc2];
                                    xp1 = (xo1+x1)*xB;
                                    xp2 = (xo2+x2)*xB;
                                    trg = [xp1(:)'; xp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,NG,NG);
                                    all = all./scl;
                                    %--
                                    check = 0;
                                    if check & lev>1
                                        temp = dir'*all*dir;
                                       % figure;imagesc(real(temp));axis square;colorbar;
                                    end

                                    all = dirlev{lev}'*all*dirlev{lev};%interpolation from N to NN

                                    if check & lev>1
                                       err = temp - all(1:2^(lev-1):end,1:2^(lev-1):end);
                                       if norm(err)>1e-4
                                           err
                                           pause;
                                       end
                                    end
                                    
                                    %scale
                                    t1 = [0:NT-1] + x1*NT;
                                    t2 = [0:NT-1] + x2*NT;
                                    src = [kc1; kc2];
                                    tp1 = (to1+x1)*xB;
                                    tp2 = (to2+x2)*xB;
                                    trg = [tp1(:)'; tp2(:)'];
                                    scl = fun(N,trg,src);              scl = reshape(scl,NT,NT);
                                    all = all.*scl;
                                    %fprintf(1,'%d %d %d %d %d %d\n', k1,k2,x1,x2,real(all(1,1)),imag(all(1,1)));
                                    %put
                                    u(t1+1,t2+1) = u(t1+1,t2+1) + all;
                                end
                            end
                        end
                        %x
                    end
                end
                %k
                NOW = NXT;
                clear NXT;
            end
            %----------------
        end
    end
end
%z


