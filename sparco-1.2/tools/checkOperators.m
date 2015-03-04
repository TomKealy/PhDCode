function checkOperators

% Initialize random number generators
randn('state',0); rand('state',1);

% Check basic operators
A = opBinary(103,407,0);             checkOperator(A);
A = opBinary(103,407,1);             checkOperator(A);

A = opBlur(103,407);                 checkOperator(A);

kernel = sin(linspace(0,9,351))';
A = opConvolve1d(kernel,1, 0);       checkOperator(A);
A = opConvolve1d(kernel,37,0);       checkOperator(A);
A = opConvolve1d(kernel,1, 1);       checkOperator(A);
A = opConvolve1d(kernel,41,1);       checkOperator(A);

A = opTranspose(opCurvelet2d(103,163,6,32)); 
checkOperator(A);

A = opDCT(103);                      checkOperator(A);
A = opDCT(128);                      checkOperator(A);

A = opDiag(137,pi);                  checkOperator(A);
A = opDiag(139,pi+sqrt(-2));         checkOperator(A);
A = opDiag(randn(141,1));            checkOperator(A);
A = opDiag(sqrt(randn(143,1)));      checkOperator(A);

A = opDirac(25);                     checkOperator(A);

A = opFFT(176);                      checkOperator(A);
A = opFFT(256);                      checkOperator(A);

A = opFFT2C(32,64);                  checkOperator(A);
A = opFFT2C(31,178);                 checkOperator(A);

A = opFFT2d(32,64);                  checkOperator(A);
A = opFFT2d(41,138);                 checkOperator(A);

A = opGaussian(55,79,0);             checkOperator(A);
A = opGaussian(56,78,1);             checkOperator(A);
A = opGaussian(57,77,2);             checkOperator(A);
A = opGaussian(58,76,3);             checkOperator(A);
A = opGaussian(59,75,4);             checkOperator(A);

A = opHaar(64);                      checkOperator(A);
A = opHaar(16,2);                    checkOperator(A);

A = opHaar2d(32,64);                 checkOperator(A);
A = opHaar2d(16,32,2);               checkOperator(A);

A = opHadamard(64);                  checkOperator(A);
A = opHadamard(58);                  checkOperator(A);

A = opHeaviside(67,0);               checkOperator(A);
A = opHeaviside(76,1);               checkOperator(A);

A = opMask(rand(127,1) > 0.3);       checkOperator(A);

A = opMatrix(randn(34,67));          checkOperator(A);
A = opMatrix(sqrt(randn(76,43)));    checkOperator(A);

p = randperm(137);
A = opRestriction(137,p(1:57));      checkOperator(A);

A = opSign(55,79,0);                 checkOperator(A);
A = opSign(56,78,1);                 checkOperator(A);
A = opSign(57,77,2);                 checkOperator(A);
A = opSign(58,76,3);                 checkOperator(A);

A = opToepGauss(55,79,'toeplitz',0); checkOperator(A);
A = opToepGauss(56,78,'toeplitz',1); checkOperator(A);
A = opToepGauss(57,77,'circular',0); checkOperator(A);
A = opToepGauss(58,76,'circular',1); checkOperator(A);

A = opToepSign(55,79,'toeplitz',0);  checkOperator(A);
A = opToepSign(56,78,'toeplitz',1);  checkOperator(A);
A = opToepSign(57,77,'circular',0);  checkOperator(A);
A = opToepSign(58,76,'circular',1);  checkOperator(A);

% Skip the following operators
% - opPadding
% - opReal
% - opSplitComplex
% - opWavelet

% Check meta operators
B1 = opSign(55,79,0);
B2 = opToepSign(57,77,'circular',0);
B3 = opFFT(176);
A  = opBlockDiag(B1,B2,B3);          checkOperator(A);

p  = randperm(101);
A  = opColumnRestrict(54,101,p(1:38),'discard');
checkOperator(A);
A  = opColumnRestrict(54,101,p(1:38),'zero');
checkOperator(A);

B  = opFFT(136);
A  = opCrop(B,36,50,0,0);            checkOperator(A);
A  = opCrop(B,136,136,0,0,1);        checkOperator(A);
A  = opCrop(B,36,50,10,10);          checkOperator(A);
A  = opCrop(B,36,50,-10,-10);        checkOperator(A);
A  = opCrop(B,36,50,10,-10);         checkOperator(A);
A  = opCrop(B,36,50,-60,-60);        checkOperator(A);
A  = opCrop(B,36,50,-60, 10);        checkOperator(A);
A  = opCrop(B,36,50,10,-60);         checkOperator(A);
A  = opCrop(B,36,50,110,110);        checkOperator(A);
A  = opCrop(B,36,50,110,10);         checkOperator(A);
A  = opCrop(B,36,50,10,110);         checkOperator(A);
A  = opCrop(B,36,50,210,210);        checkOperator(A);
A  = opCrop(B,36,50,210,10);         checkOperator(A);
A  = opCrop(B,36,50,10,210);         checkOperator(A);

B1 = opSign(55,79,0);
B2 = opToepSign(55,77,'circular',0);
B3 = opFFT(55);
A  = opDictionary(B1,B2,B3);         checkOperator(A);

B1 = opSign(55,79,0);
B2 = opToepSign(79,77,'circular',0);
B3 = opFFT(77);
A  = opFoG(B1,B2,B3);                checkOperator(A);

B1 = opSign(6,9,1);
B2 = opFFT(31);
A = opKron(B1,B2);                   checkOperator(A);

A = opTranspose(A);                  checkOperator(A);

B1 = opSign(55,79,0);
B2 = opToepSign(57,79,'circular',0);
B3 = opFFT(79);
A  = opStack(B1,B2,B3);              checkOperator(A);

B1 = opSign(64,64,0);
B2 = opToepSign(64,64,'circular',0);
B3 = opFFT(64);
A  = opSum(B1,B2,B3);                checkOperator(A);

B  = opToepSign(57,79,'circular',0);
A  = opWindowedOp(128,B,randn(57,1),20);
checkOperator(A);



function checkOperator(A)

cws = get(0,'CommandWindowSize');
str = opToString(A);

width = max(0,cws(1)-length(str)-17);
dots  = repmat('. ', 1,floor(width/2));

fprintf('Checking %s %*s', str, width, dots);

status = dottest(A,100,'quiet');

if status ~= 0
   fprintf('Failed\n');
else
   fprintf('OK\n');
end
