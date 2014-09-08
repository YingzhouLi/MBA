function grid = bfio_grid(EPS)
  NG = EPS;
  
  %grid = [0:EPS-1]/(EPS-1);
  
  grid = (cos([NG-1:-1:0]/(NG-1)*pi)+1)/2;
  
  %tmp = cos([NG-1/2:-1:1/2]/NG*pi);  tmp = tmp./max(abs(tmp));  grid = (tmp+1)/2;
  
  %tmp = (lglnodes(EPS-1)+1)/2;  grid = tmp(end:-1:1)';
  
