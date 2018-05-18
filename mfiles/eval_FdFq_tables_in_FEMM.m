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

function [F_map,SOLUTION] = eval_FdFq_tables_in_FEMM(geo,per,id_vett,iq_vett,eval_type,filemot,mat,dataIn)

nsim = round(geo.nsim_singt/2);
delta_sim = geo.delta_sim_singt;

temp1 = id_vett;
temp2 = iq_vett;
[Id,Iq] = meshgrid(temp1,temp2);

% results = [];

I = Id + 1i * Iq;
%io = abs(I);
iAmp = abs(I);
gamma = angle(I) * 180/pi;
[i0,~]=calc_io(geo,per);



% geo_ref=geo;
% per_temp=cell(5,5);
% for rr=1:size(io,1)
%     for cc=1:size(io,2)
%         per.io=io(rr,cc);
%         per.gamma=gamma(rr,cc);
%         per_temp{rr,cc}=per;
%     end
% end

currentDir=pwd(); %OCT

filemot = strrep(filemot,'.mat','.fem');

%io_femm=io*geo.Nbob;
iAmpCoil=iAmp*geo.Nbob*geo.n3phase;   %AS

%Fd = zeros(size(Id)); Fq = Fd; T = Fd; dT = Fd;
if dataIn.LossEvaluationCheck == 0
    for rr = 1:size(iAmp,1)
        parfor cc = 1:size(iAmp,2)
            [thisfilepath,dirName]=createTempDir();
            disp(['Evaluation of position I:',num2str(iAmp(rr,cc)),' gamma:',num2str(gamma(rr,cc))]);
            if isoctave()
                dest_dir=strcat(dirName,'mot0.fem');
                copyfile(filemot,dest_dir);
            else
                copyfile(filemot,[dirName,'mot0.fem']);
            end
            perTmp=per;
            perTmp.gamma=gamma(rr,cc);
            perTmp.overload=iAmp(rr,cc)/i0;
            [SOLUTION{rr,cc}] = simulate_xdeg(geo,perTmp,eval_type,dirName,'mot0.fem');
            cd(currentDir)
        end
    end
else
    for rr = 1:size(iAmp,1)
        parfor cc = 1:size(iAmp,2)
            [thisfilepath,dirName]=createTempDir();
            disp(['Evaluation of position I:',num2str(iAmp(rr,cc)),' gamma:',num2str(gamma(rr,cc))]);
            copyfile(filemot,[currentDir,'\tmp\',dirName,'\mot0.fem']);
            [SOLUTION{rr,cc}] = simulate_xdeg_IronLoss(geo,iAmpCoil(rr,cc),per.BrPP,gamma(rr,cc),eval_type,mat,per);
            cd(currentDir)
        end
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
        
        if isfield(SOL,'F')
            Fr(rr,cc) = mean(real(SOL.F));
            Ft(rr,cc) = mean(imag(SOL.F));
        end
    end
end

F_map.Id   = Id;
F_map.Iq   = Iq;
F_map.Fd   = Fd;
F_map.Fq   = Fq;
F_map.T    = T;
F_map.dT   = dT;
F_map.dTpp = dTpp;

if exist('Pfes_h')
    F_map.Pfe = Pfes_h+Pfes_c+Pfer_h+Pfer_c;
end


% save sim_mot_temp F_map