function [Nage,Catch]=runmodelF(N0,T,nyrs,F,M, ModelPar)
%RUNMODELF--run cod model using specified F
%
%[Nage,Catch]=runmodelF(N0,T,nyrs,F,M, ModelPar)
%
%Inputs
%       N0 = 1-by-9 array of initial population age structure
%       T = temperature anomaly (assumed fixed)
%     nyrs= number of years to simulate
%       F = fishing mortality
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
Nage=nans(nyrs+1,9);
Catch=Nage;
Nage(1,:)=N0(:)';
for j=1:nyrs;
    SSB = Nage(j,:)*ModelPar.weightSSB(:);
    SSB=SSB/1e3;%kg to tons
    
    %recuitment model--see recruitmentmodelcompare
    R=ModelPar.Rmodel.mfunc(SSB,T,ModelPar.Rmodel.coefs);
    
    R=R*1e3;%1,000 of fish to fish
    
    
    
    Nage(j+1,1)=R;%age 1 fish
    for k=2:4;%unfished
        Nage(j+1,k)=Nage(j,k-1)*exp(-M);
        Catch(j,k-1)=0;
    end
    for k=5:8;%unfished
        if(k==5);
            Mact = ModelPar.M4coefs(1)*T+ModelPar.M4coefs(2);
            Mact = Mact+M;
        else
            Mact=M;
        end
        Nage(j+1,k)=Nage(j,k-1)*exp(-(Mact+F));
        Catch(j,k-1)=(Nage(j,k-1)-Nage(j+1,k))*F/(F+Mact);%percent of deaths due to fishing
    end
    %age 9
    N8=Nage(j,8);
    N9=N8*exp(-(M+F));%number of new 9s
    C8=(N8-N9)*F/(F+M);
    Catch(j,8)=C8;
    N9p=Nage(j,9);%plus group
    N10p=N9p*exp(-(M+F));%number of new 9s
    C9p=(N9p-N10p)*F/(F+M);
    Catch(j,9)=C9p;
    Nage(j+1,9)=N9+N10p;
end
Catch=Catch(1:nyrs,:);