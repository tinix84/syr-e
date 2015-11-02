%% VERSIONE 20 11 2011
% FluxElab
% PostProcess Flux_UVW --> Flux_dq

% phase flux
% time sequence is POSITIVE
% physical sequence is POSITIVE

% coil positive is towards airgap 

% eliminates time = 0
Fluxa = Flux_1sim(2:end,1);
Fluxb = Flux_1sim(2:end,2);
Fluxc = Flux_1sim(2:end,3);
% rotor mechanical position (RunS_SetParCase_new)
% t = 0, fase corrente in base a gamma, posizione rotore 0
pos_0 = 0;
pos_Mov = time/1000 * Cas.n * pi/30 + pos_0;  % rotor position (rad mec)
pos_Mov = pos_Mov(2:end)';

%% posizione del rotore (rad elt)
if isfield(Mac,'th0')
    % casi successivi feb 2010 (+90 perchè Waveform in Magnet\coil è un sin)
    theta = Mac.th0 * pi/180 + pos_Mov * Mac.p;
else
    % casi precedenti feb 2010
    % d_axis Vs alpha_axis mechanical position
    theta = pi/2 + pos_Mov * Mac.p;  % mec rad
end
% dq flux
FluxAlpha = 2/3 * (Fluxa - 0.5*Fluxb - 0.5*Fluxc);
Fluxbeta  = 2/3 * sqrt(3)/2 * (Fluxb - Fluxc);

Fluxd = 2/3*( (Fluxa-0.5*Fluxb-0.5*Fluxc).*cos(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* sin(theta));
Fluxq = 2/3*( -(Fluxa-0.5*Fluxb-0.5*Fluxc).*sin(theta) + sqrt(3)/2 * (Fluxb - Fluxc) .* cos(theta));

Fluxd_0 = mean(Fluxd);
Fluxq_0 = mean(Fluxq);

