function res = fun0(N,x,k)
  nx = size(x,2);
  nk = size(k,2);
  %x = x/N; %LEXING: IMPORTANT
  %------------
  xk = (x(1,:)'*k(1,:) + x(2,:)'*k(2,:));
  kr = sqrt(k(1,:).^2 + k(2,:).^2);
  sx = (2 + sin(2*pi*x(1,:)).*sin(2*pi*x(2,:)))/3;
  tmp = (2*pi)* (xk + sx(:)*kr);
  %tmp = (2*pi)*(xk);
  %tmp = (2*pi) * (sx(:)*kr);
  res = complex(cos(tmp),sin(tmp));
  
  
  
  
  %res = cos(tmp) + i*sin(tmp);
  %res = exp(i*tmp);
