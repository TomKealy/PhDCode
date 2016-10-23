function y = sense(x,mode,opts)
%% Implementation of several sensing architectures
% 
% Usage: y = sense(x,mode,opts)
% Inputs: x - input signal
%         mode - 0 to compute y = Phi x
%                1 to compute y = Phi^T x
%         opts - struct specifying several required parameters
%           opts.sampler - options: 'I', 'matrix', 'rdmdemod', 'polyphase', 'rdmsample'
%           opts.matrix - matrix for 'matrix'
%           opts.pnseq - pn sequnce for 'rdmdemod'/'polyphase'
%           opts.subsamp - subsampling factor for 'rdmdemod'/'polyphase' (N/M) 
%           opts.idx - subsampling indices for 'rdmsample'     
%
% Outputs: y - output signal (Phi x or Phi^T x depending on mode)
%
% Most recent change - 9/16/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

%% Initialize parameters
sampler = opts.sampler;

switch sampler
    %% No subsampling
    case 'I' 
        y = x;
    
    %% Defined matrix
    case 'matrix'
        switch mode
            case 0 % Regular mode: y = Phi x
                y = opts.matrix*x;
            case 1 % Adjoint mode: y = Phi^T x
                y = ctranspose(opts.matrix)*x;
        end
    
    %% Random demodulator
    case 'rdmdemod'
        subsamp = opts.subsamp;
        pnseq = opts.pnseq;
                
        switch mode
            case 0 % Regular mode: y = Phi x
                
                N = length(x);
                M = N/subsamp;
        
                Mlarge = N - M*floor(subsamp);
                Msmall = M - Mlarge;
                Msmall1 = floor(Msmall/2);
                Msmall2 = Msmall - Msmall1;
        
                x = x.*pnseq;
                
                break1 = Msmall1*floor(subsamp);
                break2 = break1+Mlarge*ceil(subsamp);

                x1 = x(1:break1);
                x2 = x(break1+1:break2);
                x3 = x(break2+1:N);
                
                y1 = transpose(sum(reshape(x1,floor(subsamp),Msmall1),1));
                y2 = transpose(sum(reshape(x2,ceil(subsamp),Mlarge),1));
                y3 = transpose(sum(reshape(x3,floor(subsamp),Msmall2),1));
                y = [y1; y2; y3];
                
            case 1 % Adjoint mode: y = Phi^T x
                
                M = length(x);
                N = M*subsamp;
        
                Mlarge = N - M*floor(subsamp);
                Msmall = M - Mlarge;
                Msmall1 = floor(Msmall/2);
                Msmall2 = Msmall - Msmall1;
                
                break1 = Msmall1*floor(subsamp);
                break2 = break1+Mlarge*ceil(subsamp);
                
                y1 = pnseq(1:break1).*kron(x(1:Msmall1),ones(floor(subsamp),1));
                y2 = pnseq(break1+1:break2).*kron(x(Msmall1+1:Msmall1+Mlarge),ones(ceil(subsamp),1));
                y3 = pnseq(break2+1:N).*kron(x(Msmall1+Mlarge+1:M),ones(floor(subsamp),1));
                y = [y1; y2; y3];
        end
    
    %% Compressive multiplexor
    case 'polyphase'
        subsamp = opts.subsamp;
        pnseq = opts.pnseq;
        M = length(x)/subsamp;
        switch mode
            case 0 % Regular mode: y = Phi x
                y = sum(reshape(x.*pnseq,M,subsamp),2);
            case 1 % Adjoint mode: y = Phi^T x
                y = pnseq.*repmat(x,[subsamp,1]);
        end
        
    %% Random subsampling
    case 'rdmsample'
        subsamp = opts.subsamp;
        idx = opts.idx;
        switch mode
            case 0 % Regular mode: y = Phi x
                y = x(idx);
            case 1 % Adjoint mode: y = Phi^T x
                y = zeros(subsamp*length(x),1);
                y(idx) = x;
        end
end