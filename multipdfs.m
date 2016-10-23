
x = [-100:1:100];

noisepdf = normpdf(x, 0, 1);

sig1pdf = normpdf(x, 10, 1);

sig2pdf = normpdf(x, 10, 10);

sig3pdf = normpdf(x, 10, 100);

sig4pdf = normpdf(x, 10, 1000);

figure
plot(x, noisepdf, x, sig1pdf)

figure
plot(x, noisepdf, x, sig2pdf)

figure
plot(x, noisepdf, x, sig3pdf)

figure
plot(x, noisepdf, x, sig4pdf)
