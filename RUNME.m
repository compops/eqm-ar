%----------------------------------------------------
%
% Evaluation of EM, EqM and LS on AR(p) models With missing data
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

%% Initialize settings
clear
T=1250;                         % No. data points
svar.missingdata=1;             % Should data be missing?
svar.missingtype=0;             % 0 random, 1 random sequences
svar.outliers=0;                % Should outliers be present?
svar.nmax=15;                   % Maximum filter orders for A(q) and B(q)
svar.sigmae=1;                  % WN-variance for noise e
svar.outliervar=4;              % Variance of outliers
svar.order=[];                  % If empty random orders otherwise fixed
opt.miter=50;                   % No. maximum EM iterations
opt.coefdiff=1e-5;              % Minimum change of parameter vector
opt.minlldiff=1e-6;             % Minimum decrease of the log-likelihood
opt.coefdiffEqM=1e-5;		% Minmium difference in EqM coef.
svar.r = 1;
N=50;				% No. systems to generate
rates=0:0.05:0.7;		% The fractions of missing data to test

%% Main loop
j=1; 
for mm=rates
    for i=1:N
        % Generate data and save some information
        svar.rate=mm; dataOK=1; 
        while (dataOK)
            [sys(j),data(j)]=rndARmodel(svar,T); 
            if (max(data(j).y) > 8); dataOK=0; end
        end
        
        opt.na=sys(j).n(1);
        opt.initialtheta=-BuildPhi(data(j).ye,sys(j).n)\data(j).ye';
        
        % Calculate EM, EqM and standard solutions
        tic; outputEM(j)=ARemsub(sys(j),data(j),svar,opt); timeEM(j)=toc;
        tic; outputEqM(j)=AReqmsub(sys(j),data(j),svar,opt); timeEqM(j)=toc;
        tic; outputEqM2(j)=AReqmsubsimple(sys(j),data(j),svar,opt); timeEqM(j)=toc;
        tic; outputSTD(j)=ARstdsub(sys(j),data(j),svar,opt); timeSTD(j)=toc;
        j=j+1
    end
end

%% Data processing
j=1;
for mm=rates
    perMF(j).rate=mm; perMSE(j).rate=mm; perTIME(j).rate=mm;
    tmp=[outputEM([sys.rate]==mm).mfEM]; 
        perMF(j).Mem=mean(tmp); perMF(j).Vem=var(tmp);
    tmp=[outputEqM([sys.rate]==mm).mfEqM]; 
        perMF(j).Meqm=mean(tmp); perMF(j).Veqm=var(tmp);
    tmp=[outputEqM2([sys.rate]==mm).mfEqM]; 
        perMF(j).Meqm2=mean(tmp); perMF(j).Veqm2=var(tmp);
    tmp=[outputSTD([sys.rate]==mm).mfSTD];
        perMF(j).Mstd=mean(tmp); perMF(j).Vstd=var(tmp);
    tmp=sqrt([outputEM([sys.rate]==mm).mseEM]./[sys([sys.rate]==mm).n].^2); 
        perMSE(j).Mem=mean(tmp); perMSE(j).Vem=var(tmp);
    tmp=sqrt([outputEqM([sys.rate]==mm).mseEqM]./[sys([sys.rate]==mm).n].^2);
        perMSE(j).Meqm=mean(tmp); perMSE(j).Veqm=var(tmp);
    tmp=sqrt([outputEqM2([sys.rate]==mm).mseEqM]./[sys([sys.rate]==mm).n].^2);
        perMSE(j).Meqm2=mean(tmp); perMSE(j).Veqm2=var(tmp);
    tmp=sqrt([outputSTD([sys.rate]==mm).mseSTD]./[sys([sys.rate]==mm).n].^2);
        perMSE(j).Mstd=mean(tmp); perMSE(j).Vstd=var(tmp);
    tmp=timeEM([sys.rate]==mm); perTIME(j).Mem=mean(tmp);
    tmp=timeEqM([sys.rate]==mm); perTIME(j).Meqm=mean(tmp);
    tmp=timeEqM2([sys.rate]==mm); perTIME(j).Meqm2=mean(tmp);
    tmp=timeSTD([sys.rate]==mm); perTIME(j).Mstd=mean(tmp);
    j=j+1;
end

if svar.missingdata==1; save(['missingdata' N 'runs' num2str(now) '.mat']); end
if svar.missingdata==0; save(['outliers' N 'runs' num2str(now) '.mat']); end


%% Plotting
figure(223);
plotcolor.red=200; plotcolor.green=200; plotcolor.blue=240;
plotcolor.rgb=[plotcolor.red plotcolor.green plotcolor.blue]/255;

subplot(1,3,1);
%ciplot([perMF.Mem]-1.96*sqrt([perMF.Vem]/N),[perMF.Mem]+1.96*sqrt([perMF.Vem]/N),rates,plotcolor.rgb)
hold on
%ciplot([perMF.Meqm]-1.96*sqrt([perMF.Veqm]/N),[perMF.Meqm]+1.96*sqrt([perMF.Veqm]/N),rates,plotcolor.rgb)
%ciplot([perMF.Mstd]-1.96*sqrt([perMF.Vstd]/N),[perMF.Mstd]+1.96*sqrt([perMF.Vstd]/N),rates,plotcolor.rgb)
f=plot(rates,[perMF.Mem],'black',rates,[perMF.Meqm],'black-*',rates,[perMF.Meqm2],'black-+',...
    rates,[perMF.Mstd],'black--','LineWidth',2);
hold off
legend(f,'EM','EqM','EqM2','Std');

ylabel('Model fit'); xlabel('Fraction of outliers present');
if (svar.missingdata ==1); xlabel('Fraction of data missing');  end

subplot(1,3,2);
%ciplot([perMSE.Mem]-1.96*sqrt([perMSE.Vem]/N),[perMSE.Mem]+1.96*sqrt([perMSE.Vem]/N),rates,plotcolor.rgb)
hold on
%ciplot([perMSE.Meqm]-1.96*sqrt([perMSE.Veqm]/N),[perMSE.Meqm]+1.96*sqrt([perMSE.Veqm]/N),rates,plotcolor.rgb)
%ciplot([perMSE.Mstd]-1.96*sqrt([perMSE.Vstd]/N),[perMSE.Mstd]+1.96*sqrt([perMSE.Vstd]/N),rates,plotcolor.rgb)
plot(rates,[perMSE.Mem],'black',rates,[perMSE.Meqm],'black-*',...
    rates,[perMSE.Mstd],'black--','LineWidth',2);
hold off
ylabel('rmse'); xlabel('fraction of outliers present');
if (svar.missingdata ==1); xlabel('fraction of data missing');  end

subplot(1,3,3);
loglog(rates,[perTIME.Mem],'black',rates,[perTIME.Meqm],'black-*',rates,[perTIME.Meqm2],'black-+',...
    rates,[perTIME.Mstd],'black--','LineWidth',2);
ylabel('Computational time/system'); xlabel('Fraction of outliers present');
if (svar.missingdata ==1); xlabel('fraction of data missing');  end
