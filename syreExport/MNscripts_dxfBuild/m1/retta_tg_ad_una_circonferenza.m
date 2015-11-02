%retta tangente ad una circonferenza di centro alpha e beta ed una
%retta passante per Po=(x0,y0) e coeff_angolare=45°.
%retta nella forma:
function [a b c]=retta_tg_ad_una_circonferenza(alpha,beta,x0,y0)
a=x0-alpha;
b=y0-beta;
c=-(x0*(x0-alpha)+y0*(y0-beta));