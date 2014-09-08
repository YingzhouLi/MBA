function res = fun0(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:)) / N;
  tmp = (2*pi) * xk;  res = cos(tmp) + i*sin(tmp);
  
