clear all

L=512; %number of sub-bands
L0 = (L-1)/2;
m=200;%50;
positions = randi(L,[1,10]);%generate random spikes for signal

x0 = zeros(1,L); %Tx PSD
x0(positions) = 1;

figure
plot(x0)

S = randn(m,L);

Eb_N0_db = -5:15;
[r,c] = size(Eb_N0_db);
runs = 100;
mse_spg = zeros(c,runs);
for i=1:c
    for j=1:runs
        sigma = 10^(-(Eb_N0_db(i)/20))/m
        b = S*x0' + sigma*randn(1,m)';
           
               
        opts = spgSetParms('verbosity', 0);
        solution = spg_bpdn(S, b, 0.01, opts);
        %solution1 = spg_lasso(S,b, 1, opts);
        
        mse_spg(i,j) = norm(solution' - x0)/norm(x0);
        %mse_laaso(i,j) = norm(solution1' - x0)/norm(x0,2);
    end
end

 mse = sum(mse_spg,2)/runs
 figure
 plot(Eb_N0_db, mse,'-xb')
 xlabel('Eb-N0-dB')
 ylabel('mse')