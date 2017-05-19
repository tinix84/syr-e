% Copyright 2015
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power and torque versus speed curves at given Imax and Vmax 

% input:
%           fdfq_idiq_n256.mat
%           ReadParameters.m
%           ktMax_idiq.mat - kt max locus (MTPA)
%           kvMax_idiq.mat - kv max locus (MTPV)

% output:
%			torque and power versus speed profiles, given converter Imax and Vmax
%			Plim.mat and figures saved into a new subfolder

% Induction motor includes w_slip

close all, clear all, addpath mfiles, load LastPath

% load the two motor data: mat file and ReadParameters.m
[FILENAME, pathname, FILTERINDEX] = uigetfile([pathname '/*_n*.mat'], 'LOAD DATA');
load([pathname FILENAME]);
if exist([pathname 'ReadParameters.m'])
    run([pathname 'ReadParameters']);
else
    if ~exist('motor_name')
        prompt={'Motor Name','Axes type (SR or PM)','Number of pole pairs','Peak phase voltage [V]','Peak phase current [A]','Max speed [rpm]'};
        name='Operating Limits';
        numlines=1;
        defaultanswer={'MotName','SR','2','310','15','8000'};
        
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        motor_name = setup{1};
        axes_type = setup{2};
        p = eval(setup{3});
        Vmax = eval(setup{4});
        Imax = eval(setup{5});
        Nmax = eval(setup{6});
        
        save([pathname FILENAME],'p','Vmax','Imax','motor_name','axes_type','Nmax','-append');
    else
        prompt={'Peak phase voltage [V]','Peak phase current [A]','Max speed [rpm]'};
        name='Operating Limits';
        numlines=1;
        defaultanswer={num2str(Vmax),num2str(Imax),num2str(Nmax)};
        
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        Vmax = eval(setup{1});
        Imax = eval(setup{2});
        Nmax = eval(setup{3});
        
        save([pathname FILENAME],'Vmax','Imax','Nmax','-append');
    end
end

save LastPath pathname

rad2rpm=30/pi/p;

% pathname for output
pathname1 = [pathname 'Imax = ' num2str(Imax) '\'];
mkdir(pathname1)

% dq reference frame
if isempty(axes_type)
    if sum(Id(:,1)) < 0
        axes_type = 'PM';  % SPM style
    else
        axes_type = 'SR';  % IPM style
    end
end

n = 2;  % number of curves
Imax_vect = Imax * linspace(0.5,1,n);
Plim_all = cell(1,n);

% rotor bars temperature (IM motor)
if exist('Wslip','var')
    prompt={'Rotor temperature (IM only)'};
    name='Input for slip speed calculation';
    numlines=1;
    defaultanswer={num2str(Rr_temp)};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
else
    answer = 0;
end


for jj = 1:n
    
    Plim_all{jj} = calc_profiles_at_VmaxImax(pathname,FILENAME,Imax_vect(jj),answer);
    
    Plim = Plim_all{jj};
    F_A = Plim.F_A;
    Tmax = Plim.T_A;
    Pmax = max(Plim.P);
    
    plot_drive_limit_curves;    % flux weak limits: torque, power, voltage .. vs speed
   
end

% save mat\Plim Plim
save([pathname1 'Plim'],'Plim','Plim_all');

