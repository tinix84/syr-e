function Mass = calcMassCu(geo)
% 
% Mass = calcMassCu(geo)
% 
% Evaluate the copper mass of the stator

% SF - 07/07/2016

% Data
rhoCu=8940; % [kg/m^3]

% evaluation of lend

avv=geo.avv;
[r,c]=size(avv);
ini=0;
fin=0;
for ii=1:c
    if ini==0
        if avv(1,ii)==1 || avv(2,ii)==1
            ini=ii;
        end
    elseif fin==0
        if avv(1,ii)==-1 || avv(2,ii)==-1
            fin=ii;
        end
    end
end

if fin==0
    fin=c+1;
end

yq=fin-ini;

if yq==1  %concentrated winding
    lend=0.5*(geo.wt+pi*(geo.r+geo.g+geo.lt/2)*sin(pi/6/geo.p/geo.q));    % Gamba - A new PMASRM with nonconventional FS pole combination
else
    alpha=yq*2*pi/(6*geo.p*geo.q);
    lend=2*geo.lt+(geo.r+geo.g+geo.lt/2)*alpha;
    clear alpha
end

Mass=rhoCu*geo.Aslot/1e6*(geo.l+lend)/1000*(6*geo.p*geo.q); % [kg]

