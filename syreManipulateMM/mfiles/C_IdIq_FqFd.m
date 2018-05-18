%   IMPORTANT:
%   Id (n,m) varies along rows
%   Iq (n,m) varies along columns
%   Fd(n,m) Fq(n,m) go aggordingly

debug = 0;

[puntid, puntiq]=size(Fd);     % m' colonne

[NumD,NumQ]=size(Id);

FDmax=max(max(Fd))*0.8;
FDmin=min(min(Fd))*0.75;

FQmax=max(max(Fq))*0.95;
FQmin=min(min(Fq))*0.95;

fd=linspace(FDmin,FDmax,puntid);
fq=linspace(FQmin,FQmax,puntiq);

[FD,IQ]=meshgrid(fd,Iq(:,1));

% costruzione IQ_dato_fd_fq
%%%
% 1a) - inversione
% ID_dato_fd_iq(NumQ,puntid)

ID_dato_fd_iq = zeros(size(FD));
for n=1:1:NumQ,
    a = interp1(Fd(n,:),Id(n,:),FD(n,:),'PCHIP');
    ID_dato_fd_iq(n,:) = a;
    %    %% debug
    %    figure(100)
    %    plot(fd,ID_dato_fd_iq), grid on, pause
end

if debug
    fig3 = figure;
    surf(FD,IQ,ID_dato_fd_iq);grid on;
    xlabel('Fd [Vs]'),ylabel('iq [A]'),zlabel('id [A]')
    title('ID\_dato\_fd\_iq');
    pause(0.5)
end

% 2a) - interpolazione
% FQ_dato_fd_iq(NumQ,puntid)
FQ_dato_fd_iq = zeros(size(FD));
for n=1:1:NumQ,
    a=interp1(Id(n,:),Fq(n,:),ID_dato_fd_iq(n,:),'PCHIP');
    FQ_dato_fd_iq(n,:) = a;
end

if debug
    fig4 = figure;
    surf(FD,IQ,FQ_dato_fd_iq);grid on;
    xlabel('Fd [Vs]'),ylabel('iq [A]'),zlabel('Fq [Vs]')
    title('FQ\_dato\_fd\_iq');
    pause(0.5)
end

% 3a) - inversione
% IQ_dato_fd_fq(puntiq,puntid)
[FD,FQ]=meshgrid(fd,fq);
IQ_dato_fd_fq = zeros(size(FD));
for n=1:1:puntid,
    a=interp1(FQ_dato_fd_iq(:,n),IQ(:,n),FQ(:,n),'PCHIP');
    IQ_dato_fd_fq(:,n)=a;
end

if debug
    fig5 = figure;
    surf(FD,FQ,IQ_dato_fd_fq);
    grid on
    xlabel('Fd [Vs]'),ylabel('Fq [Vs]'),zlabel('iq [A]')
    title('IQ\_dato\_fd\_fq','FontSize',12,'FontWeight','bold');
    pause(0.5)
end

% costruzione ID_dato_fd_fq
%%%
[ID,FQ]=meshgrid(Id(1,:),fq);
% 1b) - inversione
% IQ_dato_id_fq(puntiq,NumD)
IQ_dato_id_fq = zeros(size(FD));
for n=1:1:NumD,
    a = interp1(Fq(:,n),Iq(:,n),FQ(:,n),'PCHIP');
    IQ_dato_id_fq(:,n) = a;
end
if debug
    fig6 = figure;
    surf(ID,FQ,IQ_dato_id_fq);grid on;
    xlabel('Id [A]'),ylabel('Fq [Vs]'),zlabel('iq [A]')
    title('IQ\_dato\_id\_fq');
    pause(0.5)
end

% 2b) - interpolazione
% FD_dato_fq_id(puntiq,NumD)
FD_dato_id_fq = zeros(size(FD));
for n=1:1:NumD,
    a=interp1(Iq(:,n),Fd(:,n),IQ_dato_id_fq(:,n),'PCHIP');
    FD_dato_id_fq(:,n) = a;
end

if debug
    fig7 = figure;
    surf(ID,FQ,FD_dato_id_fq);grid on;
    xlabel('Id [A]'),ylabel('Fq [Vs]'),zlabel('Fd [Vs]')
    title('FD\_dato\_id\_fq');
    pause(0.5)
end

% 3a) - inversione
% ID_dato_fd_fq(puntiq,puntid)
[FD,FQ]=meshgrid(fd,fq);
ID_dato_fd_fq = zeros(size(FD));
for n=1:1:puntiq,
    a=interp1(FD_dato_id_fq(n,:),ID(n,:),FD(n,:),'PCHIP');
    ID_dato_fd_fq(n,:)=a;
end

if debug
    fig8 = figure;
    surf(FD,FQ,ID_dato_fd_fq);
    grid on
    xlabel('Fd [Vs]'),ylabel('Fq [Vs]'),zlabel('id [A]')
    title('ID\_dato\_fd\_fq','FontSize',12,'FontWeight','bold');
    pause(0.5)
    figure(fig5)
end


% fdfq_idiqs to be evaluated (flux coord)
TF = 3/2 * p * (FD .* IQ_dato_fd_fq - FQ .* ID_dato_fd_fq);
FF = sqrt(FD.^2 + FQ.^2);
% delta = atan(FQ./FD);
deltaF = zeros(size(FD));
deltaF(FD ~= 0) = atan(FQ(FD ~= 0)./FD(FD ~= 0));
IF = sqrt(ID_dato_fd_fq.^2 + IQ_dato_fd_fq.^2);
ID_dato_fd_fq(ID_dato_fd_fq == 0) = eps;
argIF = atan(IQ_dato_fd_fq./ID_dato_fd_fq);
% fi = pi/2 + delta - argI;    % PF angle
% PF = cos(fi);

if isoctave()  %OCT
    name_file = strcat(pathname, 'IdIq_FdFq.mat');
    save ('-mat7-binary', name_file,'fd','fq','IQ_dato_fd_fq','ID_dato_fd_fq', ...
        'TF','FF','IF','deltaF','argIF');
    clear name_file
else
    save ([pathname 'IdIq_FdFq.mat'],'fd','fq','IQ_dato_fd_fq','ID_dato_fd_fq', ...
        'TF','FF','IF','deltaF','argIF');
end
