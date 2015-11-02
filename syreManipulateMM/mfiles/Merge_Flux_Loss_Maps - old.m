% 09 07 2014
% input:
% ldlq_idiq ... mag model n x n
% ldlq_idiq ... loss model m x m
% output
% ldlq_idiq_n256_mag_and_loss.mat : homogeneous maps 256 x 256

%%  
clear all, close all, clc
Utilizzo_TuttiIDati = 1;
% carica i dati delle matrici flussi correnti
load ultimo
[nomefile, PATHNAME] = uigetfile([PATHNAME '\*.mat'],'IdIqFdFq');
percorso = [PATHNAME nomefile];
load(percorso);
save('ultimo','PATHNAME');
% 13-01-2011 BB
Id_ModMgn = Id; Iq_ModMgn = Iq;
[nomefile_loss, PATHNAME] = uigetfile([PATHNAME '\*.mat'],'LOSS MAP');
if ischar(nomefile_loss)
    load([PATHNAME nomefile_loss]);
end
% 13-01-2011 BB
Id_Loss = Id; Iq_Loss = Iq;

% 13-01-2011 BB
% Controllo se sono nel caso: (mappa delle perdite piu' ampia della mappa magnetica)
% Se è così posso utilizzare al limite i dati (Id,Iq) della mappa magnetica
%if (min(min(Id_ModMgn))>min(min(Id_Loss)) || min(min(Iq_ModMgn))>min(min(Iq_Loss)) || max(max(Id_ModMgn))<max(max(Id_Loss)) || max(max(Iq_ModMgn))<max(max(Iq_Loss)))
%    Utilizzo_TuttiIDati=0;
%end
if Utilizzo_TuttiIDati==0
    id_min = max([min(min(Id_ModMgn)),min(min(Id_Loss))]); id_max = min([max(max(Id_ModMgn)),max(max(Id_Loss))]);
    iq_min = max([min(min(Iq_ModMgn)),min(min(Iq_Loss))]); iq_max = min([max(max(Iq_ModMgn)),max(max(Iq_Loss))]);
else
    id_min = min([min(min(Id_ModMgn)),min(min(Id_Loss))]); id_max = max([max(max(Id_ModMgn)),max(max(Id_Loss))]);
    iq_min = min([min(min(Iq_ModMgn)),min(min(Iq_Loss))]); iq_max = max([max(max(Iq_ModMgn)),max(max(Iq_Loss))]);
end
id_vett=linspace(id_min,id_max,256); iq_vett=linspace(iq_min,iq_max,256);
[Id,Iq]=meshgrid(id_vett,iq_vett);
% Flussi
Fd=interp2(Id_ModMgn,Iq_ModMgn,Fd,Id,Iq,'spline'); Fq=interp2(Id_ModMgn,Iq_ModMgn,Fq,Id,Iq,'spline');
% Perdite
if exist('Pe','var')
    Pe=interp2(Id_Loss,Iq_Loss,Pe,Id,Iq,'spline'); Pc=interp2(Id_Loss,Iq_Loss,Pc,Id,Iq,'spline'); Ph=interp2(Id_Loss,Iq_Loss,Ph,Id,Iq,'spline');
    if exist('Pes','var')
        Pes=interp2(Id_Loss,Iq_Loss,Pes,Id,Iq,'spline'); Pcs=interp2(Id_Loss,Iq_Loss,Pcs,Id,Iq,'spline'); Phs=interp2(Id_Loss,Iq_Loss,Phs,Id,Iq,'spline');
        Per=interp2(Id_Loss,Iq_Loss,Per,Id,Iq,'spline'); Pcr=interp2(Id_Loss,Iq_Loss,Pcr,Id,Iq,'spline'); Phr=interp2(Id_Loss,Iq_Loss,Phr,Id,Iq,'spline');
    end
end
if exist('P_BarRot','var')
    P_BarRot=interp2(Id_Loss,Iq_Loss,P_BarRot,Id,Iq,'spline');
end
if exist('Pmag','var')
    Pmag=interp2(Id_Loss,Iq_Loss,Pmag,Id,Iq,'spline');
end
if Utilizzo_TuttiIDati==1
    caselle_da_azzerare = ((round(Id)>max(max(Id_Loss)))+(Id<min(min(Id_Loss)))+(Iq>max(max(Iq_Loss)))+(Iq<min(min(Iq_Loss))))>0;
%     keyboard
    if exist('Pe','var')
        Pe(caselle_da_azzerare)=0; Pc(caselle_da_azzerare)=0; Ph(caselle_da_azzerare)=0;
        if exist('Pes','var')
            Pes(caselle_da_azzerare)=0; Pcs(caselle_da_azzerare)=0; Phs(caselle_da_azzerare)=0;
            Per(caselle_da_azzerare)=0; Pcr(caselle_da_azzerare)=0; Phr(caselle_da_azzerare)=0;
        end
    end
    if exist('P_BarRot','var')
        P_BarRot(caselle_da_azzerare)=0;
    end
    if exist('Pmag','var')
        Pmag(caselle_da_azzerare)=0;
    end
end

save([PATHNAME 'ldlq_idiq_n256_mag_and_loss.mat'],'Id','Iq','Fd','Fq','P*','coeff_Pfe','velDim') 
% % opzioni da utente
% % 10-11-2010 
% % Aggiunto il tipo_motore ('PM' oppure 'SR')
% prompt={'riavvolgimento (N''/N)','allungamento (L''/L)', ...
%         'lungh testate p.u.','n livelli coppia','n livelli velocità','tipo_motore','riduzione Pmag x segmentazione assiale'};
% name='Input';
% numlines=1;
% defaultanswer={'1','1','0','8','10','PM','1'};
% setup=inputdlg(prompt,name,numlines,defaultanswer);
% 
% Kr = eval(setup{1});
% Kl = eval(setup{2});
% Kt = eval(setup{3});
% nt = eval(setup{4});
% ns = eval(setup{5});
% tipo_motore = (setup{6});
% segmentazione_assiale = eval(setup{7});