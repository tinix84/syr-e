% post_proc.m
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

mi_loadsolution;

temp_out = mo_getcircuitproperties('fase1');
temp_out = temp_out - mo_getcircuitproperties('fase1n');
f1 = temp_out(3) * 2 * p;
temp_out = mo_getcircuitproperties('fase2');
temp_out = temp_out - mo_getcircuitproperties('fase2n');
f2 = temp_out(3) * 2 * p;
temp_out = mo_getcircuitproperties('fase3');
temp_out = temp_out - mo_getcircuitproperties('fase3n');
f3 = temp_out(3) * 2 * p;

ang0 = th_m; 
% evaluate torque
% T1 - linea vicina al rotore
x = xr + gap*1/6;
ang1 = 180/p + ang0;
[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(180/p,0.5);

T1 = mo_lineintegral(4);
T1 = T1(1) * 2 * p;
mo_clearcontour();

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

%% 22/01/2013 MG Decommentare le seguenti righe per eseguire il calcolo delle perdite:

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

%% Comentare la seguente riga e decommentare le righe precdenti per eseguire il calc delle Pfe
% Bx0=0; By0=0; Bx1=0; By1=0; Bx2=0; By2=0; Bx3=0; By3=0;

mo_close, mi_close

% flussi dq
fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);

%% Costruzione matrice delle soluzioni:
sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3]) Bx0 By0 Bx1 By1 Bx2 By2 Bx3 By3];

% sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];


