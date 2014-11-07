function res = fun6(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:));
  kr2 = k(1,:).^2 + k(2,:).^2;
  tmp = (2*pi)*(xk + 2*pi*ones(nx,1)*kr2) ;
  res = complex(cos(tmp),sin(tmp));
