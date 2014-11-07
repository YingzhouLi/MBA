function res = fun6(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:));
  x2k2 = (x(1,:).^2 + x(2,:).^2)'*ones(1,nk) + ones(nx,1)*(k(1,:).^2 + k(2,:).^2);
  tmp = (2*pi)*(sqrt(2)*xk + 0.5*x2k2) ;
  res = complex(cos(tmp),sin(tmp));
