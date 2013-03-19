%----------------------------------------------------
%
% Generating random AR(p) systems with missing data
%
%----------------------------------------------------
%
% Parameter inference in AR processes with missing data
%
% Authors: Johan Dahlin, Fredrik Lindsten, 
%          Thomas B. Schön
%
% Presented at ERNSI workshop.
% Maastricht, NL, 2012
%
%----------------------------------------------------


function [sys0, data] = rndARmodel(svar,T)

%% Preliminary calculations
na=1; nb=2;         % Initial orders
nest=floor(2*T/3);    % Partition point between estimation and test sets

% Generate model orders (strictly proper)
na = randi(svar.nmax);
if ~isempty(svar.order) na=svar.order; end

%% Generate filters A(q)
nom = genpolynomial(na);

% Export the true system parameters in sys0-struct.
sys0.n = na; a = poly(nom.roots);
sys0.a = a(2:end); theta = [sys0.a]';
sys0.roots=nom.roots;

%% Generate signal
% Generate noise
e = sqrt(svar.sigmae)*randn(1, 2*T); sys0.noisevar=svar.sigmae;

% Generate output
y = zeros(1, 2*T+na);

for t = 1:(2*T)
    phi = [-y(t:t+na-1)]; % -y(t-na) .. -y(t-1)
    phi = phi((na):-1:1); % Flip
    y(t+na) = phi*theta + e(t);
end

data.y = y(na+1+T:end); % discard the na initial zeros
data.e = e(1+T:end);
data.t = (1:T);

%% Missing data
if (svar.missingdata == 1) % && (svar.rate>0))
    sys0.rate=svar.rate; 
    
    if (svar.missingtype==0)
        sys0.type='random missing data'; 
        
        % Randomly select missing data points
        data.nmissing=floor(nest*svar.rate);
        data.tmissing=randsample(1+na:nest-na,data.nmissing);
        
    elseif (svar.missingtype==1)
        sys0.type='random sequences of missing data'; 
        
        % Randomly select sequences of missing data points
        data.nmissing=floor(nest*svar.rate); data.tmissing=[];
        
        % While there are more missing points to be places, randomly select
        % number of missing points and start time.
        Tmissing=0;
        while (Tmissing < data.nmissing)
            nseqmiss=rand*data.nmissing;
            tstart=ceil(rand*nest);
            if (~((tstart+nseqmiss)>(nest-na)) && (tstart>na+1))
                data.tmissing=[data.tmissing tstart:(tstart+nseqmiss)];
                Tmissing=Tmissing+nseqmiss;
            end
        end
        
        % truncate missing vector
        data.tmissing=data.tmissing(1:data.nmissing);
    end
    
    % Replace them with zeros/the noise realization
    data.ymissing=data.y(data.tmissing);
    data.y(data.tmissing)=0; %data.e(data.tmissing);
end

%% Outliers
if (svar.outliers == 1) % && (svar.rate>0))
    sys0.type='outliers'; sys0.rate=svar.rate; 
    
    % Randomly select missing data points
    data.noutlying=floor(T*svar.rate);
    data.toutlying=randsample(1+na:nest-na,data.noutlying);
    data.youtlying=data.y(data.toutlying);
    
    % Replace them with the outliers noise realization
    data.y(data.toutlying)=sqrt(svar.outliervar)*randn(data.noutlying,1);
end

%% Return data
% Divide data set into estimation and test sets
data.ye = data.y(1:nest);
data.ee = data.e(1:nest); data.te = data.t(1:nest);

data.yt = data.y(nest+1:end);
data.et = data.e(nest+1:end); data.tt = data.t(nest+1:end);
data.T=T; data.Te=nest;

%% Subroutine for filter pole generation
function polyn = genpolynomial(N)
% genpolynomial Generate a random polynomial with N roots.
% [polyn] = genpolynomial(N) returns a struct with the model order (n),
% the number of complex roots (k) and real roots (l), the roots are also
% returned in the struct (roots). The number of complex and real roots 
% are randomly selected.
%

polyn.normlimit=0.9;
polyn.n=N; % number of roots
polyn.k=2*floor(randi(polyn.n)/2); % number of complex roots
polyn.l=polyn.n-polyn.k; % number of real roots

% Generate real roots
correct=0;
%while (~correct)
    polyn.roots=-polyn.normlimit+rand(polyn.l,1)*2*polyn.normlimit;
%    if (min(abs(polyn.roots)) > (1-polyn.normlimit)); correct=1; end
%end
% Generate complex roots
for (j=1:polyn.k/2)
    correct=0;
    while (~correct)    
       % Generate point in upper unit square
       x1=-polyn.normlimit+rand*2*polyn.normlimit; 
       x2=rand*polyn.normlimit; 
       z=x1+x2*sqrt(-1);
       
        %Accept if within the unit circle
        if ((abs(z) < polyn.normlimit)) % && (abs(z) > (1-polyn.normlimit)))
            correct=1;
        end
    end
    
    % Store root with its complex conjugate
    polyn.roots=[polyn.roots; x1+x2*sqrt(-1); x1-x2*sqrt(-1)];
end
