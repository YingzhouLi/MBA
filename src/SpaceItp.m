function uOut = SpaceItp(N,m,uIn)
% 2D interpolation using FFT2
tmp = fft2(uIn);
tmp = [tmp(1:end/2,1:end/2) zeros(m/2,N-m) tmp(1:end/2,end/2+1:end);
    zeros(N-m,N);
    tmp(end/2+1:end,1:end/2) zeros(m/2,N-m) tmp(end/2+1:end,end/2+1:end)];
uOut = ifft2(tmp)*(N/m)^2;
if 0
    fprintf('check interpolation\n');
    tmp = uIn(1:5,1:5);
    tmp2 = uOut(1:2:10,1:2:10);
    [tmp(:).';tmp2(:).']
end