% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


% identificazione.m
% - simula n punti di lavoro con corrente di ampiezza io per valutare il
% modello magnetico della macchina

function [F_map,SOLUTION] = eval_FdFq_tables_in_MN(geo,per,id_vett,iq_vett,eval_type, pathname,FileName,mat,dataIn)

nsim = round(geo.nsim_singt/2);
delta_sim = geo.delta_sim_singt;

temp1 = id_vett;
temp2 = iq_vett;
[Id,Iq] = meshgrid(temp1,temp2);

%% creation folder map
Idstr=num2str(max(abs(id_vett)),3); Idstr = strrep(Idstr,'.','A');
Iqstr=num2str(max(abs(iq_vett)),3); Iqstr = strrep(Iqstr,'.','A');

if isempty(strfind(Idstr, 'A'))
    Idstr = [Idstr 'A'];
end

if isempty(strfind(Iqstr, 'A'))
    Iqstr = [Iqstr 'A'];
end

FolderName = [pathname, FileName(1:end-4) '_F_map_' Idstr 'x' Iqstr '_MN' '\'];
mkdir(FolderName);
copyfile ([pathname FileName(1:end-4) '.mn'], [FolderName]);
pathname = FolderName;  %new pathname Map
%%

I = Id + 1i * Iq;
iAmp = abs(I);
gamma = angle(I) * 180/pi;
[i0,~]=calc_io(geo,per);

currentDir=pwd(); %OCT

iAmpCoil=iAmp*geo.Nbob*geo.n3phase;   %AS

%Fd = zeros(size(Id)); Fq = Fd; T = Fd; dT = Fd;
% if dataIn.LossEvaluationCheck == 0
h = OpenMagnet(1);

for rr = 1:size(iAmp,1)
    %parfor cc = 1:size(iAmp,2)
    for cc = 1:size(iAmp,2)
        %         [thisfilepath,dirName]=createTempDir();
        disp(['Evaluation of position I:',num2str(iAmp(rr,cc)),' gamma:',num2str(gamma(rr,cc))]);
        %             if isoctave()
        %                 dest_dir=strcat(dirName,'mot0.fem');
        %                 copyfile(filemot,dest_dir);
        %             else
        %                 copyfile(filemot,[dirName,'mot0.fem']);
        
        perTmp=per;
        perTmp.gamma=gamma(rr,cc);
        perTmp.overload=iAmp(rr,cc)/i0;
        [SOLUTION{rr,cc}] = simulate_xdeg_MN(geo,perTmp,eval_type,pathname,FileName,h);
        cd(currentDir)
    end
end

Fd   = zeros(size(Id));
Fq   = zeros(size(Id));
T    = zeros(size(Id));
dT   = zeros(size(Id));
dTpp = zeros(size(Id));
if isfield(SOLUTION{1,1},'F')
    Fr = zeros(size(Id));
    Ft = zeros(size(Id));
end
if isfield(SOLUTION{1,1},'Pfes_h')
    Pfes_h = zeros(size(Id));
    Pfes_c = zeros(size(Id));
    Pfer_h = zeros(size(Id));
    Pfer_c = zeros(size(Id));
end
Ppm  = zeros(size(Id));
for rr=1:size(iAmp,1)
    for cc=1:size(iAmp,2)
        SOL=SOLUTION{rr,cc};
        dT(rr,cc) = std(SOL.T);
        dTpp(rr,cc) = max(SOL.T)-min(SOL.T);
        %ris_sim = mean(SOL(1:end,:),1);
        T(rr,cc) = -mean(SOL.T);
        Iq(rr,cc)=mean(SOL.iq);
        Id(rr,cc)=mean(SOL.id);
        Fd(rr,cc)=mean(SOL.fd);
        Fq(rr,cc)=mean(SOL.fq);
        
        if isfield(SOL,'Pfes_h')
            Pfes_h(rr,cc) = mean(SOL.Pfes_h);
            Pfes_c(rr,cc) = mean(SOL.Pfes_c);
            Pfer_h(rr,cc) = mean(SOL.Pfer_h);
            Pfer_c(rr,cc) = mean(SOL.Pfer_c);
        end
        Ppm(rr,cc) = mean(SOL.Ppm);
        if isfield(SOL,'F')
            Fr(rr,cc) = mean(real(SOL.F));
            Ft(rr,cc) = mean(imag(SOL.F));
        end
    end
end

Command='Call close(False)';
invoke(h.magnetHandler, 'processCommand', Command);


F_map.Id   = Id;
F_map.Iq   = Iq;
F_map.Fd   = Fd;
F_map.Fq   = Fq;
F_map.T    = T;
F_map.dT   = dT;
F_map.dTpp = dTpp;

if exist('Pfes_h')
    F_map.Pfe = Pfes_h+Pfes_c+Pfer_h+Pfer_c;
    F_map.Pfes_h = Pfes_h;
    F_map.Pfes_c = Pfes_c;
    F_map.Pfer_h = Pfer_h;
    F_map.Pfer_c = Pfer_c;
end

F_map.Ppm = Ppm;

F_map.speed = per.EvalSpeed;
