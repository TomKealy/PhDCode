clear all

g_sr = exprnd(1, [50,121]);
psi = zeros(1,8);
psi(2) = 1;
bet = zeros(1,121);
bet(5) = 1;

phi_r = kron(g_sr,psi)*kron(psi,bet)';
size(phi_r)