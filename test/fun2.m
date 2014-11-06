function res = fun2(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  p = bfio_k2p(N,k);
  xk = (ones(nx,1)*k(1,:)).*sqrt((p(1,:).^2)’*ones(1,nk) + (x(2,:).^2)’*(p(2,:).^2));
  tmp = (2*pi) * xk;
  res = cos(tmp) + i*sin(tmp);

