function [test_err,Ytest_learned] = elastic_net_test(Xtest,Ytest,beta,err_type)
%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%
%ELASTIC_NET_TEST evaluates the error on a test set.
%   [TEST_ERR] = ELASTIC_NET_TEST(XTEST,YTEST,BETA) returns the 
%   classification error committed by BETA on the test set XTEST, YTEST.
% 
%   [TEST_ERR,YTEST_LEARNED] = ELASTIC_NET_TEST(XTEST,YTEST,BETA) also 
%   returns the predicted values for the test samples.
%
%   [...] = ELASTIC_NET_TEST(XTEST,YTEST,BETA,ERR_TYPE) if ERR_TYPE='CLASS'
%   (default) evaluates the classification error; if ERR_TYPE='REGR'
%   evaluates the regression error.

if nargin<3, error('too few input'), end
if nargin<4; err_type = 'class'; end
if nargin>4, error('too many input'), end

Ytest_learned = Xtest*beta; % predicted values 
if isequal(err_type,'class') % evaluates classifcation error
    test_err = length(find((sign(Ytest_learned).*sign(Ytest))~=1))/length(Ytest);
elseif isequal(err_type,'regr') % evaluates regression error
    test_err = norm(Ytest_learned-Ytest)^2/length(Ytest);
end