% post_proc_v1.m

%% modifica 08 gennaio 2010 - gp
% - scarico Bx By in vari punti per valutare Pfe
% - aggiungo i vettori Bx By al fondo del vettore sol
% - Bx0, By0 sono le induzioni delle guide sull'asse di simmetrica (vicino
% ai pont radiali)
% - Bx1, By1 e Bx2,By2 sono le induzione delle guide vicino alle punte
% delle barriere

%% modifica 10 gennaio 2010 - gp
% - scarico Bx By anche nel dente di statore
% - organizzazione dei vettori x_fe, y_fe
% x_fe(1:nlay) : barriere di rotore, vicino alle punte tonde
% x_fe(nlay+1:nlay+2) : espansione dente statore (speculari)
% x_fe(nlay+3:nlay+4) : dente statore (speculari)
% x_fe(nlay+5) : giogo statore

r_fe = geo.r_fe;
hf = geo.hf;
beta_f = geo.beta_f;
l = geo.l;
p = geo.p;
x_fe = geo.x_fe;
y_fe = geo.y_fe;
Br =(4e-7*pi)*geo.Hc;        

mi_loadsolution;
% Plot flux density and save it:
 mo_showdensityplot(1,0,0,3,'mag');
 Br_string=num2str(Br,2); Br_string(Br_string=='.')='p';
 mo_savebitmap([pathname,filemot(1:end-4),'_io_',num2str(io,4),'_gamma_',num2str(gamma-90),'_Br_',num2str(Br,2) '.bmp'])
mo_hidedensityplot;
 mi_saveas([pathname,filemot(1:end-4),'_io_',num2str(io,4),'_gamma_',num2str(gamma-90),'_Br_',num2str(Br,2) '.fem'])

 % evaluate flux
temp_out = mo_getcircuitproperties('fase1');
temp_out = temp_out - mo_getcircuitproperties('fase1n');
f1 = temp_out(3) * 2 * p;
temp_out = mo_getcircuitproperties('fase2');
temp_out = temp_out - mo_getcircuitproperties('fase2n');
f2 = temp_out(3) * 2 * p;
temp_out = mo_getcircuitproperties('fase3');
temp_out = temp_out - mo_getcircuitproperties('fase3n');
f3 = temp_out(3) * 2 * p;

% evaluate torque
% T1 - linea vicina al rotore
x = xr + gap*1/6;
ang0 = th_m; ang1 = 180/p + th_m;
[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(180/p,0.5);

%mo_makeplot(2,200);
mo_makeplot(2,200,'C:\Users\matteo\Desktop\tesi\MOGA_matlab_mod\MOGA_PROT_ott11 - PARFOR - CUP 2 - Copia\Bn.txt',1);
fid = fopen('Bn.txt', 'r');
Bn = fscanf(fid, '%g %g', [2 inf]);    % It has two rows now.
Bn = Bn';
fclose(fid);
angle=Bn(:,1)/x1*180/pi*p;
B=Bn(:,2);
%figure, plot(angle,B), grid,
%
T1 = mo_lineintegral(4);
T1 = T1(1) * 2 * p;
mo_clearcontour();

% flux density on rotor magnets:
Bn_rot_e=[];
Bn_rot_i=[];
Bn_rot=[];
for ij=1:nlay
[x32r,y32r]=rot_point(geo.X32(ij),-geo.Y32(ij),pi/2/p);
mo_selectpoint(x32r,y32r);
[x32r,y32r]=rot_point(geo.X32(ij),geo.Y32(ij),pi/2/p);
mo_selectpoint(x32r,y32r);
mo_makeplot(1,200,'C:\Users\matteo\Desktop\tesi\MOGA_matlab_mod\MOGA_PROT_ott11 - PARFOR - CUP 2 - Copia\Bn.txt',1)
fid = fopen('Bn.txt', 'r');
Bn = fscanf(fid, '%g %g', [2 inf]);    % It has two rows now.
%Bn_rot_e = [Bn_rot_e,Bn'];
Bn_rot_e=Bn';
fclose(fid);
mo_clearcontour();
% keyboard
% [x42r,y42r]=rot_point(geo.X42(ij),-geo.Y42(ij),pi/2/p);
% mo_selectpoint(x42r,y42r);
% [x42r,y42r]=rot_point(geo.X42(ij),geo.Y42(ij),pi/2/p);
% mo_selectpoint(x42r,y42r);
% mo_makeplot(2,200,'C:\Users\matteo\Desktop\tesi\MOGA_matlab_mod\MOGA_PROT_ott11 - PARFOR - CUP 2 - Copia\Bn.txt',1)
% fid = fopen('Bn.txt', 'r');
% Bn = fscanf(fid, '%g %g', [2 inf]);    % It has two rows now.
% %Bn_rot_i = [Bn_rot_i,Bn'];
% Bn_rot_i=Bn';
% fclose(fid);
% mo_clearcontour();
%Bn_rot=[Bn_rot,[[mean([Bn_rot_e(:,1),Bn_rot_i(:,1)],2)],[mean([Bn_rot_e(:,2),Bn_rot_i(:,2)],2)]]];
Bn_rot=[Bn_rot,Bn_rot_e];
end
%
% T2 - linea vicina allo statore
x = xr + gap*5/6;
ang0 = -pc; ang1 = 180/p-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(180/p,0.5);
T2 = mo_lineintegral(4);
T2 = T2(1) * 2 * p;
mo_clearcontour();

% T3 - linea mediana
x = xr + gap*1/2;
ang0 = -pc; ang1 = 180/p-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(180/p,0.5);
T3 = mo_lineintegral(4);
T3 = T3(1) * 2 * p;
mo_clearcontour();

%% PERDITE NEL FERRO
% Pfe ROT - valuto il campo nelle guide di flusso
% - la guida più superficiale (numero 1) non viene considerata
% - la guida più interna (vicino all'albero) è trattata in modo particolare (punto 4 in disegna_rotore.m) perchè
% non ha spessore uniforme (nella zona dell'albero il campo è molto debole e non significativo) 
%% old geometry decommentare le righe di seguito riportate (Matteo 28/11/2011)
% Punto 0: centro guida sull'asse di simmetria del polo
%x_fe0 = geo.x0 - r_fe;
%ang_temp = (90/p + th_m)*pi/180;
%[x_fe0,y_fe0] = rot_point(x_fe0,0,ang_temp);

%% VARIAZIONE nuova geometria (Matteo 28/11/2011):
x_fe0=r_fe;
ang_temp = (90/p + th_m)*pi/180;
[x_fe0,y_fe0] = rot_point(x_fe0,0,ang_temp);

%%
% Bxy rotore
temp_Bx0 = zeros(1,nlay); temp_By0 = temp_Bx0;
temp_Bx1 = temp_Bx0; temp_By1 = temp_Bx0;
temp_Bx2 = temp_Bx0; temp_By2 = temp_Bx0;
% Bxy statore
temp_Bx3 = zeros(1,ns/2); temp_By3 = temp_Bx3;

% Punto 1: centro guida vicino al traferro
[x_fe1,y_fe1] = rot_point(x_fe(1:nlay),y_fe(1:nlay),ang_temp);

% Punto 2: specchio di punto 1
[x_fe2,y_fe2] = rot_point(x_fe(1:nlay),-y_fe(1:nlay),ang_temp);
%keyboard
% scarico il campo di rotore
for ii = 1:nlay
    temp = mo_getpointvalues(x_fe0(ii),y_fe0(ii));
    temp_Bx0(ii) = temp(2);
    temp_By0(ii) = temp(3);
    temp = mo_getpointvalues(x_fe1(ii),y_fe1(ii));
    temp_Bx1(ii) = temp(2);
    temp_By1(ii) = temp(3);
    temp = mo_getpointvalues(x_fe2(ii),y_fe2(ii));
    temp_Bx2(ii) = temp(2);
    temp_By2(ii) = temp(3);
end
[Bx0,By0] = rot_point(temp_Bx0,temp_By0,-ang_temp);
[Bx1,By1] = rot_point(temp_Bx1,temp_By1,-ang_temp);
[Bx2,By2] = rot_point(temp_Bx2,temp_By2,-ang_temp);

% Pfe STAT - valuto il campo nel primo dente

% valuto tutti i denti di statore (ns/2)
% x_fe3 = zeros(1,ns/2); y_fe3 = x_fe3;
x_fe3 = x_fe(nlay+3);
y_fe3 = y_fe(nlay+3);

for ii = 1:ns/2
    ang_cava = pc * pi/180 * (ii-1);
    [temp_x,temp_y] = rot_point(x_fe3,y_fe3,ang_cava);
    temp = mo_getpointvalues(temp_x,temp_y);
    temp_Bx3(ii) = temp(2);
    temp_By3(ii) = temp(3);
end
Bx3 = temp_Bx3;
By3 = temp_By3;

mo_close, mi_close

% flussi dq
fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);

sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3]) Bx0 By0 Bx1 By1 Bx2 By2 Bx3 By3];
mag=[angle,B];
% sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];


