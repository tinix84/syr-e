function plot_conductor_prop(mat)

BackColor=[[0.6 0.8 1.0];[0.941 0.941 0.941]];
h0=figure();
set(h0, 'Units', 'centimeters', 'Position', [2 2 15 18]);
xpos=[0.02,0.3];
ypos=0.86;
w = 0.25;
h = 0.026;

% title
uicontrol('Parent',h0, ...
        'Units','normalized', ...
        'FontSize',12, ...
        'FontName','TimesNewRoman', ...
        'FontWeight','Bold', ...
        'Style','text',...
        'Position', [0.02 0.9 0.90 0.04],...
        'BackgroundColor','w', ...
        'String',mat.MatName);
ypos = ypos-0.03;
% category
ypos = ypos-0.03;
uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(1) ypos w h], ...
          'BackgroundColor',BackColor(1,:), ...
          'String','Category:');

uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(2) ypos w h], ...
          'BackgroundColor',BackColor(2,:), ...
          'String','Conductor');

% mass density
ypos = ypos-0.03;
uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(1) ypos w h], ...
          'BackgroundColor',BackColor(1,:), ...
          'String','Mass density [kg/m^3]');

uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(2) ypos w h], ...
          'BackgroundColor',BackColor(2,:), ...
          'String',num2str(mat.kgm3));


% conductivity
ypos = ypos-0.03;
uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(1) ypos w h], ...
          'BackgroundColor',BackColor(1,:), ...
          'String','Conductivity [S/m]');

uicontrol('Parent',h0, ...
          'Units','normalized', ...
          'FontSize',9, ...
          'FontName','TimesNewRoman', ...
          'Style','text', ...
          'Position', [xpos(2) ypos w h], ...
          'BackgroundColor',BackColor(2,:), ...
          'String',num2str(mat.sigma));

