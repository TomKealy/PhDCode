function [xqhn,P,a]=prdctr(en,xqs,M,P,a,alpha)
% One-step-ahead predictor based on RLSE
% Input:  en   = Prediction error 
%         xqs  = Predictor input history [1 xq(n-M) ... xq(n-1)] 
%         M    = Order of the predictor
%         P    = Matrix to be initialized as, say, an identity matrix
%         a    = Previous parameter estimate of dimension M+1
%         alpha= Forgetting factor <=1
% Output: xqhn = Predictor output 
%         P    = Updated P matrix
%         a    = Updated parameter estimate of dimension M+1
K= P*xqs/(xqs'*P*xqs +alpha); % Predictor gain
xqhn= a'*xqs;  % Predictor output
a= a +K*en;    % Parameter update
P= (eye(M+1,M+1)-K*xqs')*P/alpha; % Parameter update
