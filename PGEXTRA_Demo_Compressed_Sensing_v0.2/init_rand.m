function seed=init_rand(seed)
%seed=10000;
if ~exist('seed','var')
    seed=floor(10000000*(now*10-floor(now*10)));
end;
tempS=RandStream('mcg16807','Seed',seed);
RandStream.setGlobalStream(tempS);%setDefaultStream(tempS);
end