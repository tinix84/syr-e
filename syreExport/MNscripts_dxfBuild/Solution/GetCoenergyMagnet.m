function cE = GetCoenergyMagnet(magData)

% GetCoenergyMagnet.m [v1.00.00 (30-11-2012)]
% Gets co-energy value
% =========================================================================
% Syntax: cE = getCoenergyMagnet(h)
% Input:
%         - magData: Magnet's data structure
%
% Output:
%         - cE:      co-energy value
% =========================================================================

mh = magData.magnetHandler;
invoke(mh,'processCommand','CALL setVariant(0,getDocument().getSolution().getCoenergy(1))');
cE = invoke(mh, 'getVariant', 0);


