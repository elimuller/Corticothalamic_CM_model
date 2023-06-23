function [nrg] = PdistGaussKern(dat,ds,bandwidth)

if nargin < 3
    bandwidth = 1;
end

%dat is msd values
pd = fitdist(dat,'Kernel','BandWidth',bandwidth);
yLC = pdf(pd,ds);
nrg = -1.*log(yLC);

end