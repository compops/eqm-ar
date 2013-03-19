%----------------------------------------------------
%
% EqM for AR(p) models With missing data
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


%% Initialize
function output = AReqmsub(sys,data,svar,opt)

Ne=length(data.ye); Nt=length(data.yt); r=0;
thetahat(1,:)=opt.initialtheta;

for k=1:opt.miter
% Eq-step
if (svar.missingdata==1); tc=data.tmissing; end
if (svar.outliers==1); tc=data.toutlying; end;
yex=[]; te=1:data.Te; te(tc)=[]; yex=data.ye(te); 

ys=data.ye; A=zeros(Ne-sys.n,Ne);
for l=1:Ne-sys.n; A(l,l:l+sys.n)=[1 thetahat(k,:) ]; end;

A=sparse(A); A1=A(1:sys.n,1:sys.n); A2=A(end-sys.n+1:end,end-sys.n+1:end).';
Rinv_matrix=A'*A;
Rinv_matrix(end-sys.n+1:end,end-sys.n+1:end) =...
    Rinv_matrix(end-sys.n+1:end,end-sys.n+1:end)+A1*A1'-A2*A2';

tmp=Rinv_matrix(tc,:); tmp2=1/det(tmp(:,tc));
if r<tmp2; r=tmp2*exp(svar.r); end
uc_k=-1*((tmp(:,tc))\(tmp(:,te)*yex'));

if ~isempty(tc)
  u=ones(length(tc),1);
  ys(1,tc)=uc_k+sqrt(log(r/tmp2)/(u'*tmp(:,tc)*u))*u; 
end

% M-step
yst=ys';
PSI=zeros(sys.n,Ne-sys.n);
for l=1:sys.n; PSI(l,:)=yst(sys.n-l+1:Ne-l).'; end;
a_est=-1*(PSI'\yst(sys.n+1:Ne));
thetahat(k+1,:)=a_est';

% Terminate?
output.breakreason='maxiter';
if ((k~=1) && (norm(thetahat(k+1,:)-thetahat(k,:)) < opt.coefdiffEqM))
    output.breakreason='coefdiff';
    break
end

end

% Model testingopt
output.initialtheta=opt.initialtheta;
output.thetahatEqM=thetahat;
output.yhatEqM=-BuildPhi(data.yt,data.ye,sys.n)*output.thetahatEqM(k+1,:)';
output.mfEqM=1-norm(data.yt'-output.yhatEqM)/norm(data.yt'-mean(data.yt'));
output.mseEqM=mean((output.thetahatEqM(k,:)-sys.a).^2);
output.y=ys;

%-----------------------------------------------
% End of File
%-----------------------------------------------

