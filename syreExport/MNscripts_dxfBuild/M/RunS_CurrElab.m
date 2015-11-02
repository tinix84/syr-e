%% VERSIONE 20 11 2011
% CurrElab
% PostPrecessed Curr (1 --> Ncoil) --> dq Curr

Ia = Curr_1sim(2:end,1);
Ib = Curr_1sim(2:end,2);
Ic = Curr_1sim(2:end,3);

% % dq Curr
% Ialpha = 2/3 * (Ia-0.5*Ib-0.5*Ic);
% Ibeta  = 2/3 * sqrt(3)/2 * (Ib - Ic);

% Currd =2/3*( (CurrU-0.5*CurrV-0.5*CurrW).*cos(theta*Mac.p*pi/180) + sqrt(3)/2 * (CurrV - CurrW) .* sin(theta*Mac.p*pi/180));
% Currq =2/3*(-(CurrU-0.5*CurrV-0.5*CurrW).*sin(theta*Mac.p*pi/180) + sqrt(3)/2 * (CurrV - CurrW) .* cos(theta*Mac.p*pi/180));


