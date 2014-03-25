%% eval_x3 - 2014 01 05

% re-evaluates the Pareto front output by the MOOA

%% pivot_cost (1 or 2, default 1)
% when re-ordering the Pareto front, sorts the soultions according to one
% of the objective functions
% pivot_cost = 1, COST_sorted is ordered according to the values of T
% pivot_cost = 2, COST_sorted is ordered according to the values of dT

clear all; close all; clc
addpath mfiles
addpath('c:\femm42\mfiles\');

if matlabpool('Size')>0
    matlabpool close force
end
matlabpool

pivot_cost = 1;

[filename, pathname_ini] = uigetfile('results\end*.mat', 'Pick a file');
load([pathname_ini filename]);
dir_name = strrep(filename,'end_','');
dir_name = strrep(dir_name,'.mat','');
pathname = [pathname_ini '\' dir_name];
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(pathname);

eval_type = 'singt';

if isempty(MESSAGE)
    runcase = 'Y';
else
    runcase = questdlg('Overwrite existing files??','Warning','Y','N','N');
end
if ~strcmp(runcase,'Y')
    return
end
[data0_name, data0_pathname] = uigetfile('data0*.m', 'Attach the appropriate data0.m');

[STATUS,MESSAGE,MESSAGEID] = copyfile([data0_pathname data0_name],'data0.m');

% if run_fem
COST = [];
% end

if exist('Pareto_front')
    x = Pareto_front;
end
%% Da eliminare
% x(:,1)=x(:,1)*mean([6,10])/mean([15,27]);
% run([data0_pathname,data0_name]);
%%
HC = zeros(size(x,1),geo.nlay); ALPHA = HC; PONT = HC; NR = zeros(size(x,1),geo.nlay+1); HF = zeros(size(x,1),geo.nlay+1); Delta_X=[];
T = zeros(size(x,1),1); dT = T; PFES = T; PFER = T; FD90 = T;
keyboard
for m = 1:size(x,1)
    
    mot_index = num2str(m);
    
    if m<10
        mot_index = ['0' mot_index];
    end
    
    % FEA evaluates each machine of the front
    cost = FEMMfitness(x(m,:)');
    
    COST = [COST; cost];
    
    % saves the bitmap of each machine
    openfemm
    opendocument('mot0.fem'); mi_savebitmap([pathname '\mot_'  mot_index '.bmp'])
    closefemm
    
    %%
    
    movefile('mot0.fem',[pathname '\mot_'  mot_index '.fem']);
    
    load geo_mot_temp
    load sim_mot_temp
    
    % geo
    save([pathname '\mot_' mot_index '.mat'],'geo','cost','out','BLKLABELS','per','rotore2','statore');
    
    % geometry
    HC(m,:) = geo.hc;
    HF(m,:) = geo.hf;
    ALPHA(m,:) = geo.alpha;
    PONT(m,:) = geo.pont;
    NR(m,:) = geo.nr;
    
    if (strcmp(geo.RotType,'3U'))
        Delta_X(m,:)=geo.Delta_X;
    end
    
    % cost functions
    T(m) = out.Tn;
    dT(m) = out.ripple_pu;
    
end

[STATUS,MESSAGE,MESSAGEID] = copyfile('data0.m',[pathname '\data0.m']);

if matlabpool('Size')>0
    matlabpool close force
end

geo_all = struct('ALPHA', ALPHA, 'HC', HC, 'HF', HF, 'PONT', PONT, 'NR', NR);

if (strcmp(geo.RotType,'3U'))
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

name_case = strrep(name_output,'sort_','')


%% Pareto front
figure(1)
hold on
for ii=1:n_mot
    plot(T(ii),dT(ii),'x'),
    text(T(ii),dT(ii)+0.1,num2str(ii));
end
grid on, hold off
xlabel('Torque [Nm]'),
ylabel('Ripple [pu]')
figure_title = strrep(filename,'.mat','');
title(figure_title)
saveas(gcf,[pathname 'pareto2x-' name_case]);

%% Bar charts
figure(2), %set(2,'Position',[40 50 1000 400],'Name',name_output)

bw = 0.7;   % bar width
set(gca,'xTickLabel',I)

subplot(3,2,1)
bar(dalpha,bw,'stacked'), grid,
legenda=[];
for kk=1:geo.nlay
    legenda{kk}=['\delta\alpha',num2str(kk)];
end
leg=legend(legenda);set(leg,'Orientation','Horizontal','Location','Best')
xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
name_output(name_output == '_') = '-';

figure(2), subplot(3,2,2)
bar(geo_all_sorted.HC,bw,'stack'), grid,
legenda=[];
for kk=1:geo.nlay
    legenda{kk}=['air (hc',num2str(kk),')'];
end
leg = legend(legenda);
set(leg,'Orientation','Horizontal','Location','Best')
xlim([0 n_mot+1]);
set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),

figure(2), subplot(3,2,3)
bar(geo_all_sorted.HF,bw,'stacked'), grid,
legenda=[];
for kk=1:geo.nlay+1
    legenda{kk}=['steel (hf',num2str(kk),')'];
end
leg = legend(legenda);
set(leg,'Orientation','Horizontal','Location','Best')
xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),

figure(2), subplot(3,2,4)
bar(x_sorted(:,end)), grid, legend('gamma'), xlim([0 n_mot+1]);
set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),

figure(2), subplot(3,2,5)
bar(abs(out_all_sorted.T)), grid, legend('Nm'); xlim([0 n_mot+1]);
set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),

subplot(3,2,6)
bar(abs(out_all_sorted.dT)), grid, legend('ripple pu');
xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),

saveas(gcf,[pathname 'bar-' name_case '_sort' num2str(pivot_cost)]);

% Cost-function riferimento figura 3 di pubblication 1:
[Y I]=sort(-dT);
COSTf=[abs(T(I,:)),dT(I,:)];

figure(3)
bar(COSTf,'Group'), grid, legend('Torque Nm','ripple pu','Orientation','Horizontal','Location','Best')
xlim([0 n_mot+1]); set(gca,'xTicklabel',I);set(gca,'xTick',1:1:length(I));
title('ripple in descending order')
saveas(gcf,[pathname 'NEW GEOMETRY Torque_ripple_fn-' name_case '_sort' num2str(pivot_cost)]);

[front,idx] = FastParetoEstimation(x,COST);
nmot_actual_front=find(idx); %macchine sul fronte vero

figure(10)
hold on

for ii=1:length(nmot_actual_front) %ii = 1:length(x) PERCHE' USARE LENGTH(X) ABBIAMO LA LUNGHEZZA DI T ???????????
    T_front(ii)=front(ii,end-1);dT_front(ii)=front(ii,end);
    plot(T_front(ii),dT_front(ii),'x'),
    text(T_front(ii),dT_front(ii)+0.1,num2str(nmot_actual_front(ii)));
end

[T_front_sorted,ii_tfs]=sort(T_front);plot(T_front(ii_tfs),dT_front(ii_tfs),'LineWidth',2)
set(gca,'FontName','Arial','FontSize',12)
grid on, hold off
xlabel('Torque - Nm'),
% ylabel('ripple - Nm')
ylabel('Torque ripple - %')
figure_title = ['3C nlay = ' num2str(geo.nlay) ' p = ' num2str(geo.p) ' - magnete plasto Br = ' num2str(geo.Br) ' T'];
title(figure_title)
saveas(gcf,[pathname 'pareto2x_actual-' name_case]);

%% Evaluation of the Progressive distribution of the Pareto Front
Pop=output.MatrixPop;
Fit=output.MatrixFitness;
legenda={};
color={'k' 'r' 'g' 'c','m'};
figure(100);
ii=1;
for jk=ceil(linspace(1,size(Fit,3),5))
    [front,idx] = FastParetoEstimation(Pop(:,:,jk),Fit(:,:,jk));
    nmot_actual_front=find(idx); %macchine sul fronte vero
    hold on;
    plot(Fit(nmot_actual_front,1,jk),Fit(nmot_actual_front,2,jk),'*','Color',color{ii});
    legenda{ii}=['front ',num2str(jk),' of ',num2str(size(Fit,3))];
    ii=ii+1;
    
end
hold off;grid on; xlabel('Torque [Nm]'); ylabel('ripple [%]'); title('Evolution of the Pareto Front during the optimization process'); legend(legenda);
saveas(gcf,[pathname 'pareto2x_evolution_during_optimization-' name_case]);
