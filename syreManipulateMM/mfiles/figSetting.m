function figSetting(width,height)
% 
% Modification: now don't change the root definitions (Simone)

if nargin()<2
    width = 12;
    height = 10;
end

% Definitions for plotting figures
% set(groot, 'DefaultTextInterpreter', 'Tex');
% set(groot, 'DefaultLegendInterpreter', 'Tex');
% set(groot, 'defaultAxesTickLabelInterpreter','Tex');

set(gcf,'defaultTextInterpreter','Latex');
set(gcf,'defaultLegendInterpreter','Latex');
set(gcf,'defaultAxesTickLabelInterpreter','Latex');

% set(gcf,'defaultFigureRenderer','painters');

set(gcf,'defaultAxesLineWidth',1);
set(gcf,'defaultLineLineWidth',1.5);

set(gcf,'defaultAxesGridLineStyle',':');
set(gcf,'defaultAxesYGrid','on');
set(gcf,'defaultAxesXGrid','on');
set(gcf,'defaultAxesZGrid','on');
set(gcf,'defaultAxesXColor',0*[1 1 1]);
set(gcf,'defaultAxesYColor',0*[1 1 1]);
set(gcf,'defaultAxesZColor',0*[1 1 1]);

set(gcf,'defaultAxesBox','on');

set(gcf,'defaultAxesNextPlot','add');

set(gcf,'defaultAxesFontSize',12);
set(gcf,'defaultTextFontSize',12);

set(gcf,'defaultAxesFontName','Times');
set(gcf,'defaultTextFontName','Times');

screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figPos(1)=screenPos(3)/2-width/2;
figPos(2)=screenPos(4)/2-height/2;
figPos(3)=width;
figPos(4)=height;

set(gcf,'Units','centimeters');
set(gcf,'Position',figPos);
set(gcf,'Color',[1 1 1]);
set(gcf,'PaperUnits','centimeter','PaperPosition',[0 0 width height])
