%----------------------------------------------------
%
% LS for AR(p) models With missing data
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



function output = ARstdsub(sys,data,svar,opt)
if (svar.outliers==1)
    m=armax(iddata(data.ye'),[opt.na 0],'LimitError',1.5);
    output.thetahatSTD=-m.ParameterVector;
end

if (svar.missingdata==1)
    output.thetahatSTD=BuildPhi(data.ye,sys.n)\data.ye';
end

output.yhatSTD=BuildPhi(data.yt,data.ye,sys.n)*output.thetahatSTD;
output.mfSTD=1-norm(data.yt'-output.yhatSTD)/norm(data.yt'-mean(data.yt'));
output.mseSTD=mean((output.thetahatSTD'-sys.a).^2);

end

%-----------------------------------------------
% End of File
%-----------------------------------------------
