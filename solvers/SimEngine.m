function SimEngine(NumSignals,Channel_axis,SNR_axis,filename,diaryfile,varargin)

%% Simulation results
if (not(isempty(diaryfile)))
    if (exist(diaryfile))
        system(['del ' diaryfile])
    end
    diary(diaryfile)
end
if (length(varargin)>=1)
    verbose = varargin{1};
else
    verbose = true;
end
if (length(varargin)>=2)
    q = varargin{2};
    if (q ~= 5)
        fprintf(1,' Cannot work on q other than 1,5.  Need to change code for that (decimation filters)');
    end
else
    q=1;
end
if (length(varargin)>=4)
    r = varargin{3};   % how many rows to draw in random
    BaseShift = varargin{4};  % The shift for non-random rows
    shiftset = true;
else
    shiftset = false;
end
if (length(varargin)>=5)
    sigtype = varargin{5};
else
    sigtype = 'sinc';
end
Results = zeros(NumSignals, length(Channel_axis)*q, length(SNR_axis));

fprintf(1, '~~~~~~~~~~~~~~~~~~~~~~~~~~ Starting Simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf(1, '~~~~~~~~~~~~   %s  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n',datestr(now));
fprintf(1, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');


%% Signal model

N = 6;                                          % Number of bands (when counting  a band and its conjugate version separately)
B = 50e6;                                   % Maximal width of each band
fnyq = 10e9;                            % Nyquist rate
Ei = [1 2 3];                               % Energy of the i'th band
Tnyq = 1/fnyq;
if (q==1)
    TimeWin = [0  10000*Tnyq];  % Time interval in which signal is observed
else
    TimeWin = [0  40000*Tnyq];  % Time interval in which signal is observed
end
Taui = [0.4 0.7 0.2]*max(TimeWin);        % Time offest of the i'th band

fprintf(1,'---------------------------------------------------------------------------------------------\n');
fprintf(1,'Signal model\n');
fprintf(1,'   N= %d, B=%3.2f MHz, fnyq = %3.2f GHz\n', N, B/1e6, fnyq/1e9);
fprintf(1,'   Ei = [%s]\n',num2str(Ei));
fprintf(1,'   Taui = [%s]\n',num2str(Taui));
fprintf(1,'   Signal type = %s\n',sigtype);

%% Sampling parameters

AliasingFactor = 195;
fp = fnyq/AliasingFactor;           % Aliasing rate
fs = q*fp;                                                    % Sampling rate at each channel, use fs=qfp, with odd q
m = max(Channel_axis);                                                      % Number of channels

% sign alternating  mixing
Mmin = 2*ceil( 0.5*fnyq/fp+0.5 ) - 1;   % The minimal value for M
%M = Mmin +(q-1);                                                     % Set to the minimal value
if (q==1)
    M = Mmin;                                                     % Set to the minimal value
else
    M = Mmin;
end
if (shiftset == false)
    r = m;                                                       % Number of random channels
    BaseShift = 1;                                    % Tap shift for nonrandom channels
end
SignPatterns = randsrc(m,M);     % draw a random +-1 for mixing sequences
for r_ind = (r+1):m                             % Update the nonrandom rows with appropriate shifts
    SignPatterns(r_ind,:) = circshift(SignPatterns(r_ind-r,:) , [ 0 BaseShift ] );
end

% calucations
Tp = 1/fp;
Ts = 1/fs;
L0 = floor(Mmin/2);                         % this is different from the paper...
L = 2*L0+1;
qtag = floor(q/2);

fprintf(1,'---------------------------------------------------------------------------------------------\n');
fprintf(1,'Sampling parameters\n');
fprintf(1,'   fp = %3.2f MHz, q=fs/fp=%d, m=%d, M=%d, r=%d, BaseShift=%d\n',fp/1e6,q,m,M,r,BaseShift);
fprintf(1,'   qtag = %d, Mmin=%d, L0 = %d, L=%d, Tp=%3.2f uSec, Ts=%3.2f uSec\n', qtag, Mmin, L0, L, Tp/1e-6, Ts/1e-6);

%% Continuous Representation
if (q==1)
    ResFactor = 5;
else
    ResFactor = 10;
end
TimeResolution = Tnyq/ResFactor;
t_axis = TimeWin(1)  : TimeResolution : TimeWin(2);

fprintf(1,'---------------------------------------------------------------------------------------------\n');
fprintf(1,'Continuous representation\n');
fprintf(1,'   Time window = [%3.2f , %3.2f) uSec\n',TimeWin(1)/1e-6, TimeWin(2)/1e-6 );
fprintf(1,'   Time resolution = %3.2f nSec, grid length = %d\n', TimeResolution/1e-9, length(t_axis));

%% Filters

% We use only even order filters to have "exactly" symmetric response around 0

% Analog filter (polyphase for speed)
dec = ResFactor * AliasingFactor / q;
h_analog = fir1(3500, 1/dec);
lenh = length(h_analog);
polylen = ceil(lenh / dec);
h_pad = [h_analog zeros(1,polylen*dec-lenh)];
polyh = reshape(h_pad,dec,polylen);

% digital expansion filter (also polyphase)
h_expand = fir1(200,1/q);
lenh_ex = length(h_expand);
polylen_ex = ceil(lenh_ex / q);
h_expand_pad = [h_expand zeros(1,polylen_ex*q - lenh_ex) ];
polyh_expand = reshape(h_expand_pad,q,polylen_ex);

% Draw all random variables so that we can debug easily later on
Carriers = zeros(NumSignals, N/2);
for SigInd=1:NumSignals
    fi = rand(1,N/2)*(fnyq/2-2*B) + B;               % Draw random carrier within [0, fnyq/2]
    Carriers(SigInd,:) = fi;
end

Sig_State_rand=cell(1);
Sig_State_randn=cell(1);
for SigInd=1:NumSignals
    state_rand=rand('state');state_randn=randn('state');
    Sig_State_rand{SigInd} = state_rand;
    Sig_State_randn{SigInd} = state_randn;
end

for SigInd=1:NumSignals
    fprintf(1,'---------------------------------------------------------------------------------------------\n');
    fprintf(1,'Generating signal %d\n',SigInd);
    %% Signal Generation
    %   fi = [1.17   4.1   2.05]*1e9;               % Draw random carrier within [0, fnyq/2]
    fi = Carriers(SigInd,:);
    x = zeros(size(t_axis));
    if (strmatch(sigtype,'sinc'))
        for n=1:(N/2)
            x = x+sqrt(Ei(n)) * sqrt(B)*sinc(B*(t_axis-Taui(n))) .* cos(2*pi*fi(n)*(t_axis-Taui(n)));
        end
    else
        Tsymbol = 1/B;
        SymbolIndex = floor(t_axis/Tsymbol);
        NumSymbols = max(SymbolIndex)+1;
        Symbols = randsrc(N/2,NumSymbols,1:4);  % 4 symbols/bit in qpsk
        for n=1:(N/2)
            x = x+sqrt(2*Ei(n)*B)*cos(2*pi*fi(n)*t_axis+(2*Symbols(n,SymbolIndex+1)-1)*pi/4);
        end
    end

    % Calculate original support set
    Sorig = [];
    % explain:  we take the starting edges: fi-B/2  and divide by fp. minus
    % 0.5 to shift half fp towards zero. then M0+1 is to move 0 = M0+1.
    Starts = ceil(  (fi-B/2)/fp-0.5+L0+1);
    Ends = ceil( (fi+B/2)/fp-0.5+L0+1);
    for i=1:(N/2)
        Sorig = union (Sorig,  Starts(i):Ends(i));
    end
    % Now add the negative frequencies
    Sorig = union(Sorig, L+1-Sorig);
    Sorig = sort(Sorig);

    %% Noise Generation
    rand('state',Sig_State_rand{SigInd});
    randn('state',Sig_State_randn{SigInd});
    noise = randn(size(t_axis));

    % Calculate energies
    NoiseEnergy = norm(noise)^2;
    SignalEnergy = norm(x)^2;
    CurrentSNR = SignalEnergy/NoiseEnergy;

    %% Mixing
    fprintf(1,'Mixing\n');

    MixedSigSequences = zeros(m,length(t_axis));
    for channel=1:m
        MixedSigSequences(channel,:) = MixSignal(x,t_axis,SignPatterns(channel,:),Tp);
    end

    MixedNoiseSequences = zeros(m,length(t_axis));
    for channel=1:m
        MixedNoiseSequences(channel,:) = MixSignal(noise,t_axis,SignPatterns(channel,:),Tp);
    end

    %% Analog low-pass filtering and actual sampling
    fprintf(1,'Filtering and decimation (=sampling)\n');


    SignalSampleSequences = [];
    NoiseSampleSequences = [];
    fprintf(1,'    Channel ');
    for channel = 1:m
        fprintf(1,'.');  if ( (mod(channel,5)==0) || (channel==m)) fprintf(1,'%d',channel); end
        SignalSequence =  MixedSigSequences(channel,:);
        NoiseSequence   =  MixedNoiseSequences(channel,:);

        [Samples3, time_axis3] = PolyFilterDecimate(SignalSequence,t_axis,polyh,lenh);
        [NSamples3, time_axis3] = PolyFilterDecimate(NoiseSequence,t_axis,polyh,lenh);

        SampleLength = length(time_axis3);
        Sample_t_axis = time_axis3;

        SignalSampleSequences(channel, 1:SampleLength) = Samples3;
        NoiseSampleSequences(channel, 1:SampleLength) = NSamples3;
    end

    fprintf(1,'\n---------------------------------------------------------------------------------------------\n');
    fprintf(1,'Sampling summary\n');
    fprintf(1,'   %d channels, each gives %d samples\n',m,SampleLength);
    %% Reconstruction stage
    fprintf(1,'---------------------------------------------------------------------------------------------\n');
    fprintf(1,'Recovery stage\n');

    % See the stages below

    %% Expand sequences, if needed

    if (q>1)
        fprintf(1,'Expanding sequences\n');
        DigitalSignalSamples = [];
        DigitalNoiseSamples   = [];
        fprintf(1,'    Channel ');
        for channel = 1:m
            fprintf(1,'.');  if ( (mod(channel,5)==0) || (channel==m)) fprintf(1,'%d',channel); end
            SignalSequence = SignalSampleSequences(channel, 1:SampleLength);
            NoiseSequence = NoiseSampleSequences(channel, 1:SampleLength);
            for factor = 1:q
                ShiftFactor  = (factor-1-qtag)*fp;
                TimeModSeq = exp(-j*2*pi* ShiftFactor *Sample_t_axis);
                SignalMod = SignalSequence .* TimeModSeq;
                NoiseMod = NoiseSequence .* TimeModSeq;
                [Samples4, Digital_time_axis] = PolyFilterDecimate(SignalMod,Sample_t_axis,polyh_expand,lenh_ex);
                [NSamples4, Digital_time_axis] = PolyFilterDecimate(NoiseMod,Sample_t_axis,polyh_expand,lenh_ex);
                location = (channel-1)*q+factor;
                DigitalLength = length(Digital_time_axis);
                DigitalSignalSamples(location, 1:DigitalLength) = Samples4;
                DigitalNoiseSamples(location, 1:DigitalLength) = NSamples4;
            end
        end
    else
        Digital_time_axis = Sample_t_axis;
        DigitalLength = SampleLength;
        DigitalSignalSamples = SignalSampleSequences;
        DigitalNoiseSamples   = NoiseSampleSequences;
    end

    fprintf(1,'\n---------------------------------------------------------------------------------------------\n');
    fprintf(1,'Expanding summary\n');
    fprintf(1,'   %d digital channels, each gives %d samples\n',m*q,DigitalLength);

    %% CTF block
    fprintf(1,'---------------------------------------------------------------------------------------------\n');
    fprintf(1,'Entering CTF block\n');

    % define matrices for fs=fp
    S = SignPatterns;
    DFTM = fft(eye(M));
    F = [ conj(fliplr(DFTM(:,2:(L0+1))))    DFTM(:,1:(L0+1)) ];
    theta = exp(-j*2*pi/L);
    np = 1:L0;
    nn = (-L0):1:-1;
    dn = [   (1-theta.^nn)./(j*2*pi*nn)      1/L    (1-theta.^np)./(j*2*pi*np)  ];
    D = diag(dn);
    A = conj(S*F*D);  % the conj comes from the revision

    % modify matrix definitions for q=fs/fp > 1
    if (q>1)
        %        F = [ F(:,(L-qtag+1):L)   F  F(:,1:qtag) ];
        F = [ conj(fliplr(DFTM(:,2:(L0+qtag+1))))    DFTM(:,1:(L0+qtag+1)) ];

        np = 1:(L0+qtag);
        nn = (-L0-qtag):1:-1;
        dn = [   (1-theta.^nn)./(j*2*pi*nn)      1/M    (1-theta.^np)./(j*2*pi*np)  ];
        D = diag(dn);
        Atemp = conj(S*F*D);
        A = zeros(m*q,L);
        for channel=1:m
            for factor = 1:q
                row = (channel-1)*q+factor;
                colstart = q+1-factor;
                A(row,:) =  Atemp(channel,  colstart:(colstart+L-1));
            end
        end
    end

    ExChannel_axis = min(Channel_axis):(max(Channel_axis)*q);

    if (verbose == true)
        fprintf(1,'---------------------------------------------------------------------------------------------\n');
        fprintf(1,'Support recovery for Signal %d\n',SigInd);
        fprintf(1,'---------------------------------------------------------------------------------------------\n');
        fprintf(1,'       | ');
        for mInd = 1:length(ExChannel_axis)
            if (ExChannel_axis(mInd) >=100)
                fprintf(1,'1');
            else
                fprintf(1,' ');
            end
        end
        fprintf(1,'\n       | ');
        for mInd = 1:length(ExChannel_axis)
            mbar = ExChannel_axis(mInd);
            digit10 = floor(mod(mbar,100)/10);
            if ( digit10*10 == mod(mbar,100) )
                fprintf(1,'%d',digit10);
            else
                fprintf(1,' ');
            end
        end
        fprintf(1,'\nSNR\\Ch | ',SigInd);
        for mInd = 1:length(ExChannel_axis)
            mbar = ExChannel_axis(mInd);
            digit1 = mod(mbar,10);
            fprintf(1,'%d',digit1);
        end
        fprintf(1,'\n-------|');
        for mInd = 1:length(ExChannel_axis)
            fprintf(1,'-');
        end
        fprintf(1,'----\n');
    end

    for SNRInd = 1:length(SNR_axis)
        SNR = SNR_axis(SNRInd);   % in dB
        SNR_val = 10^(SNR/10);          % not dB
        % combine signal and noise
        DigitalSamples = DigitalSignalSamples + DigitalNoiseSamples*sqrt(CurrentSNR/SNR_val);;

        if (verbose == true)
            fprintf(1,'%6d | ',SNR);
        end
        SuccessStr = '';
        for mInd = 1:length(ExChannel_axis)
            mbar = ExChannel_axis(mInd);

            % index the channels that participate the recovery
            RecInd = 1 : (mbar);
            Abar = A(RecInd,:);

            % Frame construction
            Q = DigitalSamples(RecInd,:)*DigitalSamples(RecInd,:)';

            % decompose Q to find frame V
            NumDomEigVals= FindNonZeroValues(eig(Q),5e-8);
            [V,d] = eig_r(Q,min(NumDomEigVals,2*N));
            v = V*diag(sqrt(d));

            % N iterations at most, since we force symmetry in the support...
            %            [u, RecSupp] = RunOMP_Unnormalized(v, Abar, N, 0, 0.01, true);
            %            RecSuppSorted = sort(unique(RecSupp));

            [RecSupp] = RunOMP_forMB(v, Abar,N,0);
            RecSuppSorted = sort(unique(RecSupp));


            % Decide on success
            if (is_contained(Sorig,RecSuppSorted)  && (rank(Abar(:,RecSuppSorted)) == length(RecSuppSorted)))
                Success= 1;
                SuccessStr = [SuccessStr 'v'];
            else
                Success = 0;
                SuccessStr = [SuccessStr '.'];
            end

            Results(SigInd,mInd,SNRInd) = Success;
        end
        if (verbose == true)
            fprintf(1,'%s\n',SuccessStr);
        end
    end

    if (mod(SigInd,25) == 0)
        save(filename)
    end

end

save(filename)

if (not(isempty(diaryfile)))
    diary off
end
