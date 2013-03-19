function Phi = BuildPhi(ys,ys0,n)
Ne=length(ys);
    
if (nargin==2)
    n=ys0;
    Phi = zeros(Ne,n);
    for (i = 1:n); Phi(i+1:end,i) = ys(1:end-i); end
end

if(nargin==3)
    Phi = zeros(Ne,n);
    for (i = 1:n)
        Phi(1:i,i)=ys0(end-i+1:end);
        Phi(i+1:end,i) = ys(1:end-i);
    end
end




