%% VERSIONE 20 11 2011
% Voltage Elab

% dq e.m.f.
Ed = Cas.n * Mac.p * pi/30 * (-Fluxq);
Eq = Cas.n * Mac.p * pi/30 * ( Fluxd);

Ed_g = Cas.n * Mac.p * pi/30 * ( Fluxq);
Eq_g = Cas.n * Mac.p * pi/30 * ( Fluxd);

pos_EMF = (pos_Mov(2:end) + pos_Mov(1:end-1))/2;
dT = (time(2) - time(1))/1000;
Ea_ =  (Fluxa(2:end)-Fluxa(1:end-1))/dT;
Eb_ =  (Fluxb(2:end)-Fluxb(1:end-1))/dT;
Ec_ =  (Fluxc(2:end)-Fluxc(1:end-1))/dT;

Ea = interp1(pos_EMF,Ea_,pos_Mov,'cubic');
Eb = interp1(pos_EMF,Eb_,pos_Mov,'cubic');
Ec = interp1(pos_EMF,Ec_,pos_Mov,'cubic');