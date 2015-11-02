
MN6 = actxserver('Magnet.application');
 MN6 = actxserver('Magnet.application');

[MachineFileName,pathname]=uigetfile([cd,'\*.mn'],'seleziona il file Magnet da usare');
set(MN6, 'Visible', 1);
DocSV = invoke(MN6, 'openDocument',[pathname,MachineFileName]);

Doc         = invoke(MN6, 'getDocument');
View        = invoke(Doc, 'getCurrentView');
Solution    = invoke(Doc, 'getSolution');
calculator  = invoke(MN6, 'getStackCalculator');

NTime       = invoke(Solution,'getFieldSolutionTimeInstants',1);  % numero istanti

Command     = 'TimeInstants = getDocument().getSolution().getFieldSolutionTimeInstants(1,time_instants)';
invoke(MN6, 'processCommand', Command);

time = [];
for ij = 1:(NTime)
    invoke(MN6, 'processCommand', ['Call setVariant(0, time_instants(' num2str(ij-1) '))']);
    time(ij) = invoke(MN6, 'getVariant', 0);
end

Coil=['U' 'V' 'W'];

Current_1sim=[];
Flux_1sim=[];
Emf_1sim=[];

for ij=1:(NTime)
    
    Command = 'ReDim ProblemID(1)';
    invoke(MN6, 'processCommand', Command);
    Command = 'ProblemID(0) = 1';
    invoke(MN6, 'processCommand', Command);
    Command = ['ProblemID(1) = ' num2str(time(ij) * 1.001)];
    invoke(MN6, 'processCommand', Command);

    
Flux_1sim_1ph = [];
Current_1sim_1ph = [];
Emf_1sim_1ph = [];

for jj=1:3
    Command = ['Call getDocument.getSolution.getFluxLinkageThroughCoil(ProblemID,' num2str(jj) ', magnitude, phase)'];
    invoke(MN6, 'processCommand', Command);
    invoke(MN6, 'processCommand', 'Call setVariant(0, magnitude)');
    z = invoke(MN6, 'getVariant', 0);
    Flux_1sim_1ph = [Flux_1sim_1ph z];
end

Flux_1sim = [Flux_1sim;Flux_1sim_1ph];

    
    for jj=1:3
        Command=['CALL getDocument().getSolution().getCurrentThroughCoil(ProblemID,"',Coil(jj),'", magnitude, phase)'];
        invoke(MN6, 'processCommand', Command);
        invoke(MN6, 'processCommand', 'Call setVariant(0, magnitude)');
        z = invoke(MN6, 'getVariant', 0);
        Current_1sim_1ph = [Current_1sim_1ph z];
    end
    Current_1sim = [Current_1sim;Current_1sim_1ph];


for jj=1:3
    Command = ['Call getDocument.getSolution.getVoltageAcrossCoil(ProblemID,' num2str(jj) ', magnitude, phase)'];
    invoke(MN6, 'processCommand', Command);
    invoke(MN6, 'processCommand', 'Call setVariant(0, magnitude)');
    z = invoke(MN6, 'getVariant', 0);
    Emf_1sim_1ph = [Emf_1sim_1ph z];
end
Emf_1sim = [Emf_1sim;Emf_1sim_1ph];
end


Fluxa = Flux_1sim(2:end,1);
Fluxb = Flux_1sim(2:end,2);
Fluxc = Flux_1sim(2:end,3);

theta=90*pi/180;

% dq current
Ia=Current_1sim(2:end,1);
Ib=Current_1sim(2:end,2);
Ic=Current_1sim(2:end,3);

IAlpha = 2/3 * (Ia - 0.5*Ib - 0.5*Ic);
Ibeta  = 2/3 * sqrt(3)/2 * (Ib - Ic);

Id = 2/3*( (Ia-0.5*Ib-0.5*Ic).*cos(theta) + sqrt(3)/2 * (Ib - Ic) .* sin(theta));
Iq = 2/3*( -(Ia-0.5*Ib-0.5*Ic).*sin(theta) + sqrt(3)/2 * (Ib - Ic) .* cos(theta));

% dq flux
FluxAlpha = 2/3 * (Fluxa - 0.5*Fluxb - 0.5*Fluxc);
Fluxbeta  = 2/3 * sqrt(3)/2 * (Fluxb - Fluxc);

Fluxd = 2/3*( (Fluxa-0.5*Fluxb-0.5*Fluxc).*cos(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* sin(theta));
Fluxq = 2/3*( -(Fluxa-0.5*Fluxb-0.5*Fluxc).*sin(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* cos(theta));

Fluxd_0 = mean(Fluxd);
Fluxq_0 = mean(Fluxq);

%% Voltage elaboration
Cas.n=1500;
Ed = Cas.n * Mac.p * pi/30 * (-Fluxq);
Eq = Cas.n * Mac.p * pi/30 * ( Fluxd);

Ea = Emf_1sim(2:end,1);
Eb = Emf_1sim(2:end,2);
Ec = Emf_1sim(2:end,3);

EAlpha = 2/3 * (Ea - 0.5*Eb - 0.5*Ec);
% EAlpha = Cas.n * Mac.p * pi/30 * (FluxAlpha);


figure;plot(time(2:end)',EAlpha,time(2:end)',IAlpha);

