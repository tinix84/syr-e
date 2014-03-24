
%% %%%%%%%%%%%%%%%
%% 2013/08/22 MG script che esegue le simulazioni in termini di singt e singm la funzione singp non è ancora stata implementata...
% la macchina deve già essere stata costruita nel formato.fem, viene
% INPUT: motXX.fem ;
% OUTPUT: simulazione puntuale, mappa;
%% %%%%%%%%%%%%%%%

clear all; close all; clc;
current_path=cd;
[pathstr, name, ext] = fileparts(current_path);

% addpath C:\Programmi\femm42\mfiles\;
% addpath('C:\Program Files (x86)\femm42\mfiles');
addpath('d:\femm42_beta\mfiles\');
addpath([pathstr,'\dxf_conv_fun\']);
addpath([pathstr,'\m2\']);
addpath([pathstr,'\m1\']);

%% domande
load ultimo.mat;

[filemot, pathname] = uigetfile([pathname,'\*.fem'], 'Evaluate FEMM - CHOOSE A MOTOR');
load([pathname,filemot(1:end-4),'_ParMachine.mat']);

choice=questdlg('Select the post process','choose post process operation', ...
    'Map','I<gamma>','I<gamma0-90>','I<gamma');
switch choice
    case 'Map'
        eval_type='singm';
        temp = inputdlg({'Nturns per slot';'Current must be a vector[A]?'},'INPUT',1,{'16','[0:2:10]'});
        eval(['Mac.N_cond=' temp{1} ';']);
        eval(['io_input=' temp{2} ';']);
        %         eval(['gamma_input=' temp{3} ';']);
        
    case 'I<gamma>'
        eval_type='singt';
        temp = inputdlg({'Nturns per slot';'Current scalar or vector[A]?';'Current phase scalar or vector[A] (dim equal to current)?'},'INPUT',1,{'16','[10 20]','[56 65]'});
        eval(['Mac.N_cond=' temp{1} ';']);
        eval(['io_input=' temp{2} ';']);
        eval(['gamma_input=' temp{3} ';']);
        
    case 'I<gamma0-90>'
        eval_type='singt';
        temp = inputdlg({'Nturns per slot';'Current must be a scalar[A]?'},'INPUT',1,{'16','10'});
        eval(['Mac.N_cond=' temp{1} ';']);
        eval(['io_input=' temp{2} ';']);
        eval(['gamma_input=' temp{3} ';']);
        
    otherwise
        break; error('UnKnown method');
end

io=io_input*Mac.N_cond;


%%  %%%%%%%%%%%%%%%%%%%%%%%%

if matlabpool('Size')>0
    matlabpool close force
end

% matlabpool(3);
% keyboard
%%  %%%%%%%%%%%%%%%%%%%%%%%
openfemm;
opendocument([pathname,filemot]);
mi_saveas([pathstr,'\mot0.fem']);
per.Vdc=300;
Mac.ns=Mac.Q/Mac.p;
Mac.Hc=0;
Mac.pathstr=pathstr;
Mac.pathname=pathstr;
% geo.TipoMot='Bru';

if strcmp(eval_type,'singt')
    
    Mac.nsim_singt = 31; Mac.delta_sim_singt = Mac.p*360/Mac.Q; gamma=gamma_input;
    
    for i = 1:length(io)
        out = valuta_motore(Mac,per,io(i),gamma(i));
        %% 2013/08/26 MG si riassegnano correnti e flussi per tenere conto del numero di spire.
        out.SOL(:,2)=out.SOL(:,2)/Mac.N_cond;
        out.SOL(:,3)=out.SOL(:,3)/Mac.N_cond;
        out.SOL(:,4)=out.SOL(:,4)*Mac.N_cond;
        out.SOL(:,5)=out.SOL(:,5)*Mac.N_cond;
        
        save('sim_mot_temp', 'out','-append');
        Istr=num2str(io_input(i),3);Istr(Istr=='.')='A';
        gammaStr=num2str(gamma,4); gammaStr(gammaStr=='.')='G';
        FILENAME = [filemot(1:end-4) '_Sim_',Istr,'_',gammaStr];
        [success,message,messageid] = mkdir(pathname,FILENAME);
        NewDir=[pathname,FILENAME,'\'];
        copyfile('sim_mot_temp.mat',[NewDir,FILENAME,'.mat']);
        geo=Mac;
        save([NewDir,FILENAME,'.mat'],'geo','-append');
    end
elseif strcmp(eval_type,'singm')
    
    Mac.nsim_singm = 16; Mac.delta_sim_singm = Mac.p*360/Mac.Q;
    ivect=io;
    %     ivect = linspace(0,io,11);
    identificazione1(Mac,ivect,ivect,0);
    %% 2013/08/26 MG caricare e riassegnare F_map serve per tenere conto del numero di conduttori in cava che è posto pari a 1
    load sim_mot_temp F_map
    F_map.Id = F_map.Id./Mac.N_cond;
    F_map.Iq = F_map.Iq./Mac.N_cond;
    F_map.Fd = F_map.Fd.*Mac.N_cond;
    F_map.Fq = F_map.Fq.*Mac.N_cond;
    save sim_mot_temp F_map;
    FILENAME = [filemot(1:end-4) '_Map_',num2str(min(io_input)),'_',num2str(max(io_input)),'_',num2str(length(io)),'x',num2str(length(io))];
    [success,message,messageid] = mkdir(pathname,FILENAME);
    NewDir=[pathname,FILENAME,'\'];
    copyfile('sim_mot_temp.mat',[NewDir,FILENAME,'.mat']);
    
    
    
    disp('Unkonwn method');
    return;
    
end


closefemm;
matlabpool close;
