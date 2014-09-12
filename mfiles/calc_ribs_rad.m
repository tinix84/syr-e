% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

%% 
[xrot_traf,yrot_traf]=calc_intersezione_cerchi(xr, r, x0);
Dx=(xr-xrot_traf(1))/5;  %per avere almeno 5 divisioni;
xcir_plot=[xr:-Dx:xr*cos(pi/2/p)];
ycir_plot=sqrt(xr^2-xcir_plot.^2);
VectCir=find(xcir_plot>=xrot_traf(1));
x_ext_rot=xcir_plot(VectCir);
y_ext_rot=ycir_plot(VectCir);

A=[];
for ii=1:nlay
    if ii==1
        X=[B2k(ii), XpBar2(ii), xxD2k(ii),xpont(ii),fliplr(x_ext_rot)];
        Y=[0, YpBar2(ii), yyD2k(ii),ypont(ii),fliplr(y_ext_rot)];
%         figure(100);hold on;fill(X,Y,'r');hold off;
        A(ii)=polyarea(X,Y);
        clear X Y;
    else
        X=[B1k(ii-1), XpBar1(ii-1), xxD1k(ii-1),xpont(ii-1),xrot_traf(ii-1),xrot_traf(ii),xpont(ii),xxD2k(ii),XpBar2(ii),B2k(ii)];
        Y=[0, YpBar1(ii-1), yyD1k(ii-1),ypont(ii-1),yrot_traf(ii-1),yrot_traf(ii),ypont(ii),yyD2k(ii),YpBar2(ii),0];
%         figure(100);hold on;fill(X,Y,'r'); hold off;
        A(ii)=polyarea(X,Y);
        clear X Y;
    end
end
Afe=cumsum(A);

rG=(B1k+B2k)/2;
M_Fe = 2*Afe*l * 1e-9 * 7800 ;   % massa ferro appeso ai ponticelli

F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;

pont = F_centrifuga/(sigma_max * l);    % mm

for jj=1:nlay
if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
    pont(jj)=0;
end    
end


hpont=pont;
rac_pont=abs(B1k-B2k)/4;

for ii=1:nlay
    
    if hpont(ii)>0
        XpontRadBarSx(ii)=B1k(ii);
        YpontRadBarSx(ii)=hpont(ii)+rac_pont(ii);
        XpontRadBarDx(ii)=B2k(ii);
        YpontRadBarDx(ii)=hpont(ii)+rac_pont(ii);
        XpontRadDx(ii)=B2k(ii)-rac_pont(ii);
        YpontRadDx(ii)=hpont(ii);
        XpontRadSx(ii)=B1k(ii)+rac_pont(ii);
        YpontRadSx(ii)=hpont(ii);

    else
        
        
        XpontRadBarSx(ii)=B1k(ii);
        YpontRadBarSx(ii)=0;
        XpontRadBarDx(ii)=B2k(ii);
        YpontRadBarDx(ii)=0;
        XpontRadDx(ii)=NaN;
        YpontRadDx(ii)=0;
        XpontRadSx(ii)=NaN;
        YpontRadSx(ii)=0;
    end
    
end

