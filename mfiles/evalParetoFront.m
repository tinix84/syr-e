
% re-evaluates all the machines of the Pareto front in 15 rotor
% positions (instead of the 5 positions with random offset used by MODE)

% pivot_cost: selects how to sort the machines of the front
% pivot_cost = 1 (default), COST_sorted is ordered according to the values of T
% pivot_cost = 2, COST_sorted is ordered according to the values of dT
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

function evalParetoFront(filename)

pivot_cost = 1;

syreRoot = fileparts(which('MODEstart.m'));
if nargin<1
    [filename, pathname_ini] = uigetfile('results\.mat', 'Pick a file');
else
    pathname_ini=[syreRoot filesep 'results' filesep];
    filename=[filename,'.mat'];
end
load(filename);
dir_name = strrep(filename,'end_','');
dir_name = strrep(dir_name,'.mat','');
pathname = [dir_name];
[~,MESSAGE,~] = mkdir(pathname);

eval_type = 'singt';

if isempty(MESSAGE)
    runcase = 'Yes';
else
    runcase = questdlg('Overwrite existing files??','Warning','Yes','No','No');
end
if ~strcmp(runcase,'Yes')
    return
end

geo=geo0;       % assign to geo intial geometric data (same in data0)
per_temp=per;   % re-assign because matlabpool don't work...
COST = [];

% filter duplicate solutions
[PSet,iA] = unique(OUT.PSet,'rows')
x = PSet;
%%

HC = zeros(size(x,1),geo.nlay); ALPHA = HC; PONT = HC; NR = zeros(size(x,1),geo.nlay+1); HF = zeros(size(x,1),geo.nlay+1); Delta_X=[];
T = zeros(size(x,1),1); dT = T; PFES = T; PFER = T; FD90 = T;

parfor m = 1:size(x,1)
    
    mot_index = num2str(m);
    disp(['Evaluating machine ' mot_index ' of ' num2str(size(x,1))])
         
    % FEA evaluates each machine of the front
    [cost{m},geometry{m},~,~] = FEMMfitness(x(m,:)',geo,per_temp,eval_type);

end
    
for m=1:size(x,1)         
    clear geo
    mot_index = num2str(m);
    if m<10
        mot_index = ['0' mot_index];
    end
    
    geo=geometry{m};
    if exist([syreRoot filesep 'empty_case.fem'],'file')>1
        empty_case_path = [syreRoot filesep 'empty_case.fem'];
    else
        empty_case_path = ['c:' filesep 'empty_case.fem'];      %TODO: fix this with something more robust and crossplatform
    end

    openfemm
    [geo] = draw_motor_in_FEMM(geo,eval_type);
    %% end drawing procedure....    
    % Save png file mot
    mi_zoomnatural;
    mi_savebitmap([pathname '\mot_'  mot_index '.bmp']);
    closefemm
    movefile([syreRoot,'\mot0.fem'],[pathname '\mot_'  mot_index '.fem']);
    % Save data geometry mot
    geo.RQ = x(m,:);
    save([pathname '\mot_' mot_index '.mat'],'geo','cost','per');
        %% %%%%%%%
    COST=[COST; cost{m}];
    if m<10
        mot_index = ['0' mot_index];
    end
    % geometry
    HC(m,:) = geo.hc;
%     HF(m,:) = geo.hf;
    ALPHA(m,:) = geo.alpha;
    PONT(m,:) = geo.pont0;
    NR(m,:) = geo.nr;
    
    % cost functions
    T(m,1) = COST(m,1); %out.Tn;
    dT(m,1) = COST(m,2); %out.ripple_pu;
    close;
end

[STATUS,MESSAGE,MESSAGEID] = copyfile([pathname_ini,filename],[pathname,'\']);

geo_all = struct('ALPHA', ALPHA, 'HC', HC, 'HF', HF, 'PONT', PONT, 'NR', NR);

if (strcmp(geo.RotType,'Seg'))
    geo_all = struct('ALPHA', ALPHA, 'HC', HC, 'HF', HF, 'PONT', PONT, 'NR', NR,'Delta_X',Delta_X);
end

out_all = struct('T', T, 'dT', dT); %, 'PFES', PFES, 'PFER', PFER, 'FD90', FD90);

%% orders the solutions of the Front according to one of the goals (T or dT)
n_mot = size(COST,1);
[Y,I] = sort(COST(:,pivot_cost));

COST_sorted = [COST(I,:) I];
x_sorted = x(I,:);

ALPHA_s = geo_all.ALPHA(I,:);
HC_s = geo_all.HC(I,:);
HF_s = geo_all.HF(I,:);
PONT_s = geo_all.PONT(I,:);
NR_s = geo_all.NR(I,:);
% DFE_s=geo_all.DFE(I,:);

T_s = out_all.T(I,:);
DT_s = out_all.dT(I,:);

geo_all_sorted = struct('ALPHA', ALPHA_s, 'HC', HC_s, 'HF', HF_s, 'PONT', PONT_s, 'NR', NR_s); %,'DFE',DFE_s,'ERROR_CODE',ERROR_CODE_s);
out_all_sorted = struct('T', T_s, 'dT', DT_s);

alpha = geo_all_sorted.ALPHA;
dalpha = diff(alpha,1,2);
dalpha = [alpha(:,1) dalpha];

name_output = strrep(filename,'end_','sort_');
name_output = strrep(name_output,'.mat','');

name_case = strrep(name_output,'sort_','');

%% Pareto front
close
figure(1), hold on
for ii=1:n_mot
    plot(T(ii),dT(ii),'x'),
    text(T(ii),dT(ii)+0.1,num2str(ii));
end
grid on, hold off
xlabel('Torque [Nm]'), ylabel('Ripple [pu]')

[front,idx] = FastParetoEstimation(x,COST);
nmot_actual_front=find(idx); %macchine sul fronte vero

figure(1), hold on
for ii=1:length(nmot_actual_front) %ii = 1:length(x) PERCHE' USARE LENGTH(X) ABBIAMO LA LUNGHEZZA DI T ???????????
    T_front(ii)=front(ii,end-1);dT_front(ii)=front(ii,end);
    plot(T_front(ii),dT_front(ii),'ro','LineWidth',2),
    text(T_front(ii),dT_front(ii)+0.1,num2str(nmot_actual_front(ii)));
end
[T_front_sorted,ii_tfs]=sort(T_front);
plot(T_front(ii_tfs),dT_front(ii_tfs),'r','LineWidth',2)
set(gca,'FontName','Arial','FontSize',12)
grid on, hold off
xlabel('Torque - Nm'),
ylabel('Torque ripple - %')
figure_title = [geo.RotType ' nlay = ' num2str(geo.nlay) ' p = ' num2str(geo.p) ' - magnete plasto Br = ' num2str(geo.Br) ' T'];
title(figure_title)
saveas(gcf,[pathname '\Pareto - ' name_case '.fig']);

%% Bar charts
figure(2), 
bw = 0.7;   % bar width
set(gca,'xTickLabel',I);
subplot(3,2,1)
bar(dalpha,bw,'stacked'), grid,
legenda=[];
for kk=1:geo.nlay
    legenda{kk}=['\Delta\alpha',num2str(kk)];
end
leg=legend(legenda);set(leg);
xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
ylabel('Mechanical Deg');
name_output(name_output == '_') = '-';

figure(2), subplot(3,2,3)
bar(geo_all_sorted.HC,bw,'stack'), grid,
legenda=[];
for kk=1:geo.nlay
    legenda{kk}=['air (hc',num2str(kk),')'];
end
leg = legend(legenda);
set(leg);
xlim([0 n_mot+1]);
set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I));
ylabel('mm');

figure(2), subplot(3,2,5)
bar(x_sorted(:,end)), grid, legend('Current Angle \gamma'), xlim([0 n_mot+1]);
set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),ylabel('Electrical Deg')

[Y I]=sort(-dT);
COSTf=[abs(T(I,:)),dT(I,:)];

figure(2), subplot(3,2,2:2:6)
bar(COSTf,'Group'), grid, legend('Torque Nm','ripple pu')
xlim([0 n_mot+1]); set(gca,'xTicklabel',I);set(gca,'xTick',1:1:length(I));
title('ripple in descending order')
saveas(gcf,[pathname '\BarChart - ' name_case '_sort' num2str(pivot_cost) '.fig']);

%% Evaluation of the Progressive distribution of the Pareto Front

% if OUT.eval_type=='MO_GA'
%     Pop=OUT.PSet;
%     Fit=OUT.PFront;
% else
%     Pop=OUT.MatrixPop;
%     Fit=OUT.MatrixFitness(:,:,end);
% end
% legenda={};
% color={'k' 'r' 'g' 'c','m'};
% figure(100);
% ii=1;
% for jk=ceil(linspace(1,size(Fit,3),5))
%     [front,idx] = FastParetoEstimation(Pop(:,:,jk),Fit(:,:,jk));
%     nmot_actual_front=find(idx); %macchine sul fronte vero
%     hold on;
%     plot(Fit(nmot_actual_front,1,jk),Fit(nmot_actual_front,2,jk),'*','Color',color{ii});
%     legenda{ii}=['front ',num2str(jk),' of ',num2str(size(Fit,3))];
%     ii=ii+1;
%     
% end
% hold off;grid on; xlabel('Torque [Nm]'); ylabel('ripple [%]'); title('Evolution of the Pareto Front during the optimization process'); legend(legenda);
% saveas(gcf,[pathname '\pareto2x_evolution_during_optimization-' name_case '.fig']);

if OUT.eval_type=='MO_OA'
    save(filename,'front','idx','-append');
    movefile([pathname],['results\']);
end
