function [B,Y,N0]=getyield(N0,T,weight,nyrs,F,M,ModelPar)
%GETYIELD--get yield, biomass, state at fixed point of codmodel
%
% [B,Y,N0]=getyield(N0,T,weight,nyrs,F,M,ModelPar)
%
% See getreferencepoints for a description of the inputs
%
% Andrew Pershing, Gulf of Maine Research Institute, 2015
%
Bp=0;B=1;
while(abs((B-Bp)/Bp)>1e-6)%find something like a fixed point
    [Nage,Catch]=runmodelF(N0,T,nyrs,F,M,ModelPar);
    B=Nage(end,:)*weight(:);
    Bp=Nage(end-1,:)*weight(:);
    Y=Catch(end-1,:)*weight(:);
    N0=Nage(end,:);
end