function [geo, volume_mag_pu] = calc_volumi(geo)

global xy_ferro

p = geo.p;                      % paia poli

Ar = geo.Ar;                    % raggio albero             
xr = geo.xr;                    % Raggio del rotore al traferro
l = geo.l;                      % lunghezza pacco

volume_tot = ((pi*((xr^2)-(Ar^2)))*1e-6) * (l*1e-3);

mi_analyze
mi_loadsolution

mo_selectblock(xy_ferro(1),xy_ferro(2));
volume_ferro = ((mo_blockintegral(5))*2*p)*(l*1e-3);
mo_clearblock()

volume_magneti =  volume_tot - volume_ferro;

volume_mag_pu = volume_magneti / volume_tot;

geo.VolumeFerro = volume_ferro;
geo.VolumeMagneti = volume_magneti;


