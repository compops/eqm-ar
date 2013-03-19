%----------------------------------------------------
%
% EM for AR(p) with missing data
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

function output = ARemsub(sys,data,svar,opt)

thetahat(1,:)=opt.initialtheta;

if svar.outliers==1; data.tmissing=data.toutlying; end
ys=data.ye;

for kk=1:opt.miter
    % E-step
    % Estimate missing data using maxmization of the conditional likelihood

    if ((svar.rate>0) && ((svar.missingdata==1)||(svar.outliers==1)))
        term1=1+sum(thetahat(kk,:).^2);
        
        for i=sort(data.tmissing)    
            term2=0;

            for k=1:sys.n
                tmp=thetahat(kk,k)';
                for j=1+k:sys.n; tmp=tmp-thetahat(kk,j)'*thetahat(kk,j-k)'; end
                term2=term2+tmp*(ys(i+k)+ys(i-k));
            end

            ys(i)=term2/term1;
        end
    end
    yy(kk,:)=ys;
    ll(kk)=sum(data.ye'-BuildPhi(ys,sys.n)*thetahat(kk,:)').^2;
    
    % M-step
    % Estimate theta using LS
    thetahat(kk+1,:)=BuildPhi(ys,sys.n)\ys';
    
    if ~(kk==1)
            output.breakreason='maxiter';
            
        if abs(ll(kk)-ll(kk-1)) < opt.minlldiff
            output.breakreason='lldiff';
            break
        end
        
        if norm(thetahat(kk+1,:)-thetahat(kk,:)) < opt.coefdiff
            output.breakreason='coefdiff';
            break
        end
        
    end
end


% Model testing
output.initialtheta=opt.initialtheta;
output.ll=ll;
output.thetahatEM=thetahat(kk+1,:);
output.yhatEM=BuildPhi(data.yt,data.ye,sys.n)*output.thetahatEM';
output.mfEM=1-norm(data.yt'-output.yhatEM)/norm(data.yt'-mean(data.yt'));
output.mseEM=mean((-output.thetahatEM-sys.a).^2);
output.y=ys;
end

%-----------------------------------------------
% End of File
%-----------------------------------------------
