function g_code=gray_code0(g_code)
N=length(g_code); N2=N/2;
if N>=4, N2=N/2; g_code(N2+1:N)=fftshift(g_code(N2+1:N)); end
if N>4, g_code=[gray_code0(g_code(1:N2)) gray_code0(g_code(N2+1:N))]; end