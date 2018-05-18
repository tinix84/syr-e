function DatiOpt = PuntiDiOttimo(d,filename)

load(filename); % fdfq_idiq_n256.mat

Vbus=d.Vbus;
p=d.p;
velmin=d.velmin;
%velbase=d.velbase;
velmax=d.velmax;

velDim=0;
ns=d.ns;
nt=d.nt;
motorType=d.motorType;


IDtot=Id(1,:);    % necessario per compatibilita'
IQtot=Iq(:,1)';   % necessario per compatibilita'

fdfq.Id=Id;
fdfq.Iq=Iq;
fdfq.Fd=Fd;
fdfq.Fq=Fq;

% torque [Nm]
%T=3/2*p*(Fd.*Iq-Fq.*Id);
T=abs(T);

% torque range
Tmax=d.Tmax;
Tmin=d.Tmin;

velmec = linspace(velmin,velmax,ns);
Tmap    = linspace(Tmin,Tmax,nt);  % [Nm]


n_coppie=max(size(Tmap));
n_vel=max(size(velmec));

% current vector
I=Id+1j*Iq;
% flusx linkage vector
F=Fd+1j*Fq;
% flux linkage amplitude
Fo=abs(F);

wmecc=velmec*pi/30;
%rpm2rads=pi/30;
%rads2Hz=1/(2*pi);
%rpm2Hz=rpm2rads*rads2Hz;

% variables init
IdMin=NaN*ones(n_vel,n_coppie);
IqMin=NaN*ones(n_vel,n_coppie);
FdMin=NaN*ones(n_vel,n_coppie);
FqMin=NaN*ones(n_vel,n_coppie);
VoMin=NaN*ones(n_vel,n_coppie);
IoMin=NaN*ones(n_vel,n_coppie);
CosfiMin=NaN*ones(n_vel,n_coppie);
P_min=NaN*ones(n_vel,n_coppie);
Potenza=NaN*ones(n_vel,n_coppie);
Coppia=NaN*ones(n_vel,n_coppie);

PERDITE_JOULE=NaN*ones(n_vel,n_coppie);
PERDITE_FERRO=NaN*ones(n_vel,n_coppie);
PERDITE_BARRE_ROTORE=NaN*ones(n_vel,n_coppie);

PERDITE_MECC=NaN*ones(n_vel,n_coppie);

BACK_EMF=NaN*ones(n_vel,n_coppie);
IRON_CURRENT=NaN*ones(n_vel,n_coppie);
SLIP=NaN*ones(n_vel,n_coppie);
ROTOR_CURRENT=NaN*ones(n_vel,n_coppie);
RES_STAT=NaN*ones(n_vel,n_coppie);

EffMOT=NaN*ones(n_vel,n_coppie);

figure()
if ~isoctave()
    figSetting;
end
if d.velmin~=d.velmax
    xlim([d.velmin d.velmax])
end

if d.Tmin~=d.Tmax
    ylim([d.Tmin d.Tmax])
end

hr=plot(0,0,'r.');
hg=plot(0,0,'g.');

disp('Map evaluation in progress...')
% per ogni velocità
for n = 1:n_vel
    for t=1:n_coppie
        %disp([int2str(n) 'of' int2str(n_vel) ' / ' int2str(t) 'of' int2str(n_coppie)])
        if Tmap(t)>=0 %motore
            % evaluation of the stator frequency
            if strcmp(motorType,'IM')
                Wslip=IM.Wslip*(1+0.004*(d.temp-d.temp0));
                FreqElet=(Wslip+velmec(n)*p*pi/30)/(2*pi);
            else
                FreqElet=(velmec(n)*p*pi/30)/(2*pi);
            end
            
            % evaluation of Pfe from model
            Pfe=calcIronLoss(d.IronLossModel,fdfq,FreqElet);
            
            % evaluation of copper rotor loss (if IM)
            if strcmp(motorType,'IM')
                Prot=3/2*IM.RR.*IM.Ir.^2;
                Prot=Prot*(1+0.004*(d.temp-d.temp0));
            else
                Prot=zeros(size(Id));
            end
            
            % back emf
            Vind=1j*2*pi*FreqElet.*F;
            
            % current component representing Fe and PM loss (Ife_m)
            Ife=2/3*Pfe./conj(Vind);
            %Ife_m=2*(Pfe)./(3*abs(Vind));
            %Ife=Ife_m.*exp(1j*angle(Vind));
            
            Io=I+Ife;
            
            % Motor voltage Voc = back emf + RI
            % Voc > Vbus is penalized by augmented loss, so to be excluded
            Rs=calcStatorResistance(d,FreqElet);
            Rs=Rs.*ones(size(Id));
            
            Vof=Vind+Rs.*Io;    % phase voltage
            Voc=sqrt(3)*Vof;    % line voltage
            cosfi=cos(angle(Io)-angle(Vof));    % PF

            % calcola i moduli
            Voc_m=abs(Voc);
            Io_m= abs(Io);
            
            % perdite meccaniche
            Pmech=polyval(d.pMechLoss,velmec(n));

            % total loss (3/2 goes with peak current values)
            Perdite=Pfe+Prot+3/2*Rs*d.n3phase.*Io_m.^2+Pmech;

            % (lim) equals 1 if Voc_m>Vbus and 0 otherwise
            lim=zeros(size(Id));
            %lim=(sign(Voc_m-Vbus)+1)/2;
            %lim(floor(Voc_m-Vbus)==0)=0;
            lim(Voc_m>Vbus)=1;
            lim(Io_m>d.Imax)=1;

            % augment loss if (lim)=(1)
            PerditeNew=Perdite+lim*1e50;
            
            if Tmap(t)<max(max(T))
                Curve=contourc(IDtot,IQtot',T,[Tmap(t) Tmap(t)]);
                idIso=Curve(1,2:end);
                iqIso=Curve(2,2:end);
                PerditeIso=interp2(Id,Iq,PerditeNew,idIso,iqIso);
                
                % calcola la coppia di correnti che porta alla minima perdita
                [PminIso, indice]=min(PerditeIso);
                % assegna le correnti al percorso ottimo solo se il limite
                % di tensione e' rispettato

                id=idIso(indice);
                iq=iqIso(indice);
                limIso=interp2(Id,Iq,lim,id,iq);
            else
                % mappe di coppia troppo piccole
                limIso=1;
            end
            
            if limIso<0.01
                %plot(velmec(n),Tmap(t),'g.')
                xdata=get(hg,'XData');
                ydata=get(hg,'YData');
                xdata=[xdata velmec(n)];
                ydata=[ydata Tmap(t)];
                set(hg,'XData',xdata,'YData',ydata);
                pause(0.01)
%                 ido = interp2(Id,Iq,real(Io),id,iq,'spline');
%                 iqo = interp2(Id,Iq,imag(Io),id,iq,'spline');
%                 if and(n>1,strcmp(motorType,'SR'))
%                     id = abs(ido);
%                 else
%                     id = ido;
%                 end
%                 iq = iqo;
                
                IdMin(n,t)=id;
                IqMin(n,t)=iq;
                FdMin(n,t)=interp2(Id,Iq,Fd,id,iq);
                FqMin(n,t)=interp2(Id,Iq,Fq,id,iq);
                VoMin(n,t)=interp2(Id,Iq,Voc_m,id,iq);
                IoMin(n,t)=interp2(Id,Iq,Io_m,id,iq);
                CosfiMin(n,t)=interp2(Id,Iq,cosfi,id,iq);
                P_min(n,t)=PminIso;
                Coppia(n,t)=Tmap(t);
                Potenza(n,t)=Tmap(t)*wmecc(n);
                PERDITE_JOULE(n,t)=interp2(Id,Iq,(3/2.*Rs*d.n3phase.*Io_m.^2),id,iq);
                PERDITE_FERRO(n,t)=interp2(Id,Iq,Pfe,id,iq);
                PERDITE_MECC(n,t)=Pmech;
                PERDITE_BARRE_ROTORE(n,t)=interp2(Id,Iq,Prot,id,iq);
                BACK_EMF(n,t)=interp2(Id,Iq,abs(Vind),id,iq);
                IRON_CURRENT(n,t)=interp2(Id,Iq,abs(Ife),id,iq);
                if exist('IM')
                    SLIP(n,t)=interp2(Id,Iq,Wslip./(2*pi*FreqElet),id,iq);
                    ROTOR_CURRENT(n,t)=interp2(Id,Iq,IM.Ir,id,iq);
                else
                    SLIP(n,t)=NaN;
                    ROTOR_CURRENT(n,t)=NaN;
                end
                RES_STAT(n,t)=interp2(Id,Iq,Rs,id,iq);
                EffMOT(n,t)=Potenza(n,t)/(Potenza(n,t)+P_min(n,t))*100;
            else
                %plot(velmec(n),Tmap(t),'r.')
                if ~isoctave()
                    xdata=get(hr,'XData');
                    ydata=get(hr,'YData');
                    xdata=[xdata velmec(n)];
                    ydata=[ydata Tmap(t)];
                    set(hr,'XData',xdata,'YData',ydata);
                    pause(0.01)
                end
            end
            
        elseif Tmap(t)<0 %generatore
            %%%
            % evaluation of the stator frequency
            if strcmp(motorType,'IM')
                Wslip=IM.Wslip*(1+0.004*(d.temp-d.temp0));
                FreqElet=(-Wslip+velmec(n)*p*pi/30)/(2*pi);
            else
                FreqElet=(velmec(n)*p*pi/30)/(2*pi);
            end
            
            % evaluation of Pfe from model
            Pfe=calcIronLoss(d.IronLossModel,fdfq,FreqElet);
            
            % evaluation of copper rotor loss (if IM)
            if strcmp(motorType,'IM')
                Prot=3/2*IM.RR.*IM.Ir.^2;
                Prot=Prot*(1+0.004*(d.temp-d.temp0));
            else
                Prot=zeros(size(Id));
            end
            
            % back emf
            if strcmp(motorType,'PM')
                Vind=1j*2*pi*FreqElet.*F;
            else
                Vind=1j*2*pi*FreqElet.*conj(F);
            end
            
            % current component representing Fe and PM loss (Ife_m)
            Ife=2/3*Pfe./conj(Vind);
            %Ife_m=2*(Pfe)./(3*abs(Vind));
            %Ife=Ife_m.*exp(1j*angle(Vind));
            
            if strcmp(motorType,'PM')
                Io=I+Ife;
            else
                Io=conj(I)+Ife;
            end
            
            % Motor voltage Voc = back emf + RI
            % Voc > Vbus is penalized by augmented loss, so to be excluded
            Rs=calcStatorResistance(d,FreqElet);
            Rs=Rs.*ones(size(Id));
            
            Vof=Vind+Rs.*Io;    % phase voltage
            Voc=sqrt(3)*Vof;    % line voltage
            cosfi=cos(angle(Io)-angle(Vof));    % PF

            % calcola i moduli
            Voc_m=abs(Voc);
            Io_m= abs(Io);
            
            % perdite meccaniche
            Pmech=polyval(d.pMechLoss,velmec(n));

            % total loss (3/2 goes with peak current values)
            Perdite=Pfe+Prot+3/2*Rs*d.n3phase.*Io_m.^2+Pmech;

            % (lim) equals 1 if Voc_m>Vbus and 0 otherwise
            lim=zeros(size(Id));
            %lim=(sign(Voc_m-Vbus)+1)/2;
            %lim(floor(Voc_m-Vbus)==0)=0;
            lim(Voc_m>Vbus)=1;
            lim(Io_m>d.Imax)=1;

            % augment loss if (lim)=(1)
            PerditeNew=Perdite+lim*1e50;
            
            if strcmp(motorType,'PM')
                flag=Tmap(t)>=min(min(T));
            else
                flag=abs(Tmap(t))<=max(max(T));
            end
            
            if flag
                if strcmp(motorType,'PM')
                    Curve=contourc(IDtot,IQtot',T,[Tmap(t) Tmap(t)]);
                else
                    Curve=contourc(IDtot,IQtot',T,abs([Tmap(t) Tmap(t)]));
                end
                idIso=Curve(1,2:end);
                iqIso=Curve(2,2:end);
                PerditeIso=interp2(Id,Iq,PerditeNew,idIso,iqIso);
                
                % calcola la coppia di correnti che porta alla minima perdita
                [PminIso, indice]=min(PerditeIso);
                % assegna le correnti al percorso ottimo solo se il limite
                % di tensione e' rispettato

                id=idIso(indice);
                iq=iqIso(indice);
                limIso=interp2(Id,Iq,lim,id,iq);
            else
                % mappe di coppia troppo piccole
                limIso=1;
            end
            
            if limIso<0.01
                %plot(velmec(n),Tmap(t),'g.')
                xdata=get(hg,'XData');
                ydata=get(hg,'YData');
                xdata=[xdata velmec(n)];
                ydata=[ydata Tmap(t)];
                set(hg,'XData',xdata,'YData',ydata);
                pause(0.01)
%                 ido = interp2(Id,Iq,real(Io),id,iq,'spline');
%                 iqo = interp2(Id,Iq,imag(Io),id,iq,'spline');
%                 if and(n>1,strcmp(motorType,'SR'))
%                     id = abs(ido);
%                 else
%                     id = ido;
%                 end
%                 iq = iqo;
                
                IdMin(n,t)=id;
                FdMin(n,t)=interp2(Id,Iq,Fd,id,iq);
                if strcmp(motorType,'PM')
                    IqMin(n,t)=iq;
                    FqMin(n,t)=interp2(Id,Iq,Fq,id,iq);
                else
                    IqMin(n,t)=-iq;
                    FqMin(n,t)=-interp2(Id,Iq,Fq,id,iq);
                end
                VoMin(n,t)=interp2(Id,Iq,Voc_m,id,iq);
                IoMin(n,t)=interp2(Id,Iq,Io_m,id,iq);
                CosfiMin(n,t)=interp2(Id,Iq,cosfi,id,iq);
                P_min(n,t)=PminIso;
                Coppia(n,t)=Tmap(t);
                Potenza(n,t)=Tmap(t)*wmecc(n);
                PERDITE_JOULE(n,t)=interp2(Id,Iq,(3/2.*Rs*d.n3phase.*Io_m.^2),id,iq);
                PERDITE_FERRO(n,t)=interp2(Id,Iq,Pfe,id,iq);
                PERDITE_MECC(n,t)=Pmech;
                PERDITE_BARRE_ROTORE(n,t)=interp2(Id,Iq,Prot,id,iq);
                BACK_EMF(n,t)=interp2(Id,Iq,abs(Vind),id,iq);
                IRON_CURRENT(n,t)=interp2(Id,Iq,abs(Ife),id,iq);
                if exist('IM')
                    SLIP(n,t)=interp2(Id,Iq,Wslip./(2*pi*FreqElet),id,iq);
                    ROTOR_CURRENT(n,t)=interp2(Id,Iq,IM.Ir,id,iq);
                else
                    SLIP(n,t)=NaN;
                    ROTOR_CURRENT(n,t)=NaN;
                end
                RES_STAT(n,t)=interp2(Id,Iq,Rs,id,iq);
                EffMOT(n,t)=Potenza(n,t)/(Potenza(n,t)+P_min(n,t))*100;
            else
                %plot(velmec(n),Tmap(t),'r.')
                if ~isoctave()
                    xdata=get(hr,'XData');
                    ydata=get(hr,'YData');
                    xdata=[xdata velmec(n)];
                    ydata=[ydata Tmap(t)];
                    set(hr,'XData',xdata,'YData',ydata);
                    pause(0.01)
                end
            end
            
        else % Tmap(t)==0: non incluso nella mappa solitamente
            
        end
    end
end

disp('Maps evaluated')

T_top_W=max(Coppia');
T_top_W(T_top_W<0)=0;
T_bot_W=min(Coppia');
T_bot_W(T_bot_W>0)=0;
T_top_W(isnan(T_bot_W))=0;
T_bot_W(isnan(T_bot_W))=0;



% Salvataggio in DatiOpt
DatiOpt.T_top_W  = T_top_W;
DatiOpt.T_bot_W  = T_bot_W;
DatiOpt.IdMin    = IdMin;
DatiOpt.IqMin    = IqMin;
DatiOpt.FdMin    = FdMin;
DatiOpt.FqMin    = FqMin;
DatiOpt.P_min    = P_min;
DatiOpt.CosfiMin = CosfiMin;
DatiOpt.VoMin    = VoMin;
DatiOpt.IoMin    = IoMin;
DatiOpt.Potenza  = Potenza;
DatiOpt.Coppia   = Coppia;

DatiOpt.PERDITE_JOULE        = PERDITE_JOULE;
DatiOpt.PERDITE_FERRO        = PERDITE_FERRO;
DatiOpt.PERDITE_MECC         = PERDITE_MECC;
DatiOpt.PERDITE_BARRE_ROTORE = PERDITE_BARRE_ROTORE;

DatiOpt.Tmap   = Tmap;
DatiOpt.velmec = velmec;
DatiOpt.IDtot  = IDtot;
DatiOpt.IQtot  = IQtot;

DatiOpt.BACK_EMF      = BACK_EMF;
DatiOpt.IRON_CURRENT  = IRON_CURRENT;
DatiOpt.SLIP          = SLIP;
DatiOpt.ROTOR_CURRENT = ROTOR_CURRENT;
DatiOpt.RES_STAT      = RES_STAT;

DatiOpt.EffMOT=EffMOT;
    



