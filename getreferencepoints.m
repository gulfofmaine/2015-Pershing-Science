function [Bmsy,Fmsy, MSY]=getreferencepoints(T,weight,M,ModelPar)
%GETREFERENCEPOINTS--get Bmsy, Fmsy for a cod model
%
% [Bmsy,Fmsy, MSY]=getreferencepoints(T,weight,M,ModelPar)
%
% Inputs
%       T = temperature anomaly
%  weight = length-9 array containing the weight of each age class (0 if
%           not included in SSB or the fishery)
%       M = assumed natural mortality
%  ModelPar = structure with parameters for temperature-dependent
%        recruitment model.  Must have the fields
%        .M4coefs = [slope, intercept] of temperature-mortality function
%             for age 4 mortality
%        .weightSSB = length-9 array of the weights for computing SSB (0
%                     implies that the age class does not contribute)
%        .Rmodel = structure describing the recruitment model
%
%  Rmodel.type    = description of the model (text)
%       .formula = formula used for the fitting (text)
%       .coefs   = coeficients from the fit
%       .stats   = [r2, p, AIC, N, r2_R, p_R]
%                  r2_R and p_R are the r2 and p statistics comparing
%                  modeled R and actual R during the out of sample period
%       .mfunc   = function handle R=mfunc(SSB,T,coefs) that returns R given SSB
%                  and temperatures.
% 
% Andrew Pershing, Gulf of Maine Research Institute, 2015
%
nyrs=10;
nsearch=15;
Fmax=2;
N0orig=ones(1,9);
N0=N0orig;
YFB=nans(nsearch+2,3);
YFB(1,2)=0;
YFB(2,2)=Fmax;
[YFB(1,3),YFB(1,1),N0]=getyield(N0,T,weight,nyrs,YFB(1,2),M,ModelPar);
if(YFB(1,3)==0);N0=N0orig;end
[YFB(2,3),YFB(2,1),N0]=getyield(N0,T,weight,nyrs,YFB(2,2),M,ModelPar);
if(YFB(2,3)==0);N0=N0orig;end
YFB(3,2)=mean(YFB(1:2,2));
[YFB(3,3),YFB(3,1),N0]=getyield(N0,T,weight,nyrs,YFB(3,2),M,ModelPar);
if(YFB(3,3)==0);
    N0=N0orig;
    YFB(3,2)=1e-4;%just a little bit
    [YFB(3,3),YFB(3,1),N0]=getyield(N0,T,weight,nyrs,YFB(3,2),M,ModelPar);
end
if(YFB(2,3)==0);N0=N0orig;end
YFB=sortrows(YFB,[2,1,3]);
j=4;
n=3;
while(j<nsearch*2+1)
    [y,I]=max(YFB(:,1));%find the current max yield
    %fprintf('%3d\t%6.4f\t%6.4f\n',[I,span(YFB(:,2))]);
    if(I~=1);
        F=mean(YFB(I+(-1:0),2));
    else
        F=mean(YFB(1:2),2);
    end
    YFB(j,2)=F;
    [YFB(j,3),YFB(j,1),N0]=getyield(N0,T,weight,nyrs,F,M,ModelPar);
    if(YFB(j,3)==0);N0=N0orig;end
    if(I<n);%look to the right
        j=j+1;
        F=mean(YFB(I+(0:1),2));
        YFB(j,2)=F;
        [YFB(j,3),YFB(j,1),N0]=getyield(N0,T,weight,nyrs,F,M,ModelPar);
        if(YFB(j,3)==0);N0=N0orig;end
        n=n+1;
    end
    n=n+1;
    j=j+1;
    YFB=sortrows(YFB,[2,1,3]);
end
[y,I]=max(YFB(:,1));%find the current max yield
Bmsy=YFB(I,3);
Fmsy=YFB(I,2);
MSY=YFB(I,1);