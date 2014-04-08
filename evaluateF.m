%% evaluateF.m
% J  [OUT] : The objective Vector. J is a matrix with as many rows as
%            trial vectors in X and as many columns as objectives.
% X   [IN] : Decision Variable Vector. X is a matrix with as many rows as
%            trial vector and as many columns as decision variables.
% Dat [IN] : Parameters defined in NNCparam.m
%
%% Main call
function [J]=evaluateF(X,Dat)

numberOfIndividuals = Dat.XPOP;

J=zeros(numberOfIndividuals,Dat.NOBJ);
f=Dat.CostProblem;

    parfor i=1:numberOfIndividuals
        
        J(i,:)=f(X(i,:)');
    end
end

%% DTLZ2 Benchmark function. Defined in:
% K. Deb, L. Tiele, M. Laummans, and E. Zitzler. Scalable test problems
% for evolutionary multi-objective optimization. Institut fur Technische
% Informatik und Kommunikationsnetze, ETH Zurich, Tech. Rep. TIK-Technical
% Report No. 112, Feb. 2001.
function J=DTLZ2(X,Dat)

Xpop=size(X,1);
Nvar=Dat.NVAR;
M=Dat.NOBJ;
K=Nvar+1-M;
J=ones(Xpop,M);

for xpop=1:Xpop
    Gxm=(X(xpop,M:Nvar)-0.5*ones(1,K))*(X(xpop,M:Nvar)-0.5*ones(1,K))';
    Cos=cos(X(xpop,1:M-1)*pi/2);
    
    J(xpop,1)=prod(Cos)*(1+Gxm);
    for nobj=1:M-1
        J(xpop,nobj+1)=(J(xpop,1)/prod(Cos(1,M-nobj:M-1)))...
            *sin(X(xpop,M-nobj)*pi/2);
    end
end
end