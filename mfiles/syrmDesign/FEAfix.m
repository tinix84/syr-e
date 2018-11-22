function [FEAfixOut] = FEAfix(dataSet,geo,map,FEAfixN)
% 
% [kd,kq] = FEAfix(dataSet,geo,map,setup)
% 

disp('FEAfix calibration...')
tic
debug=0;

xx=map.xx;
bb=map.bb;
[m,n]=size(xx);

[~, ~, ~, per, mat] = data0(dataSet);

switch FEAfixN
    case 1
        xRaw=mean(mean(xx));
        bRaw=mean(mean(bb));
    case 4
        xRaw=[xx(1,1),xx(1,end),xx(end,1),xx(end,end)];
        bRaw=[bb(1,1),bb(1,end),bb(end,1),bb(end,end)];
    case 5
        xRaw=[xx(1,1),xx(1,end),xx(end,1),xx(end,end),mean(mean(xx))];
        bRaw=[bb(1,1),bb(1,end),bb(end,1),bb(end,end),mean(mean(bb))];
    case 1000
        xRaw=reshape(xx,1,numel(xx));
        bRaw=reshape(bb,1,numel(bb));
    otherwise
        error('Put a correct number!!!')      
end
kdRaw=ones(size(xRaw));
kqRaw=ones(size(xRaw));
errorFlag=zeros(size(xRaw));

index=0;
for ii=1:geo.nlay
    RQnames{index+ii}='hc';
end
index=length(RQnames);
for ii=1:geo.nlay
    RQnames{index+ii}='dx';
end
index=length(RQnames);
RQnames{index+1}='r';
RQnames{index+2}='wt';
RQnames{index+3}='lt';
RQnames{index+4}='gamma';

geo.RQnames=RQnames;

RQ=zeros(length(RQnames),length(xRaw));
fdFEA=zeros(1,length(xRaw));
fqFEA=zeros(1,length(xRaw));
fdMod=zeros(1,length(xRaw));
fqMod=zeros(1,length(xRaw));


OBJnames{1} = 'Torque';
geo.OBJnames=OBJnames;
per.objs = [per.min_exp_torque 1];

geo.nsim_singt=4;
geo.delta_sim_singt=60;
per.overload=dataSet.CurrLoPP;
eval_type='singt';

geo.K_mesh=5;
geo.K_mesh_MOOA=10;

for mot=1:length(xRaw)
    disp(['Machine design ' int2str(mot) ' of ' int2str(length(xRaw))])
    geo.x=xRaw(mot);
    geo.b=bRaw(mot);

    r = geo.x*geo.R;                                        % rotor radius [mm]
    wt = interp2(xx,bb,map.wt,geo.x,geo.b);                       % tooth width
    wt = round(wt*100)/100;
    lt=interp2(xx,bb,map.lt,geo.x,geo.b);                         % slot length
    lt=round(lt*100)/100;
    geo.x0=geo.R * geo.x /cos(pi/2/geo.p);
    Ar(mot)=interp2(xx,bb,map.Ar,geo.x,geo.b);                          % shaft radius [mm]
    Ar(mot)=round(Ar(mot)*100)/100;
    geo.la=interp2(xx,bb,map.la,geo.x,geo.b);                         % total insulation
    
    % current phase angle
    temp_id = interp2(xx,bb,map.id,geo.x,geo.b);         % id [A]
    temp_iq = interp2(xx,bb,map.iq,geo.x,geo.b);         % iq [A]
    gamma   = round(atand(temp_iq/temp_id)*100)/100;
    
    if (isnan(lt)||(temp_iq==0))
        errorFlag(mot)=1;
    else
        hcTmp=zeros(m,n);
        dxTmp=zeros(m,n);
        for ii=1:geo.nlay
            for mm=1:m
                for nn=1:n
                    hcTmp(mm,nn)=map.hc_pu{mm,nn}(ii);
                    dxTmp(mm,nn)=map.dx{mm,nn}(ii);
                end
            end
            geo.hc_pu(ii)=interp2(xx,bb,hcTmp,geo.x,geo.b);
            geo.dx(ii)=interp2(xx,bb,dxTmp,geo.x,geo.b);
        end
        
        % fill RQ
        RQ(:,mot)=[geo.hc_pu geo.dx r wt lt gamma]';
        
        % model flux linkages
        fdMod(mot) = interp2(xx,bb,map.fd,geo.x,geo.b); % fd [Vs]
        fqMod(mot) = interp2(xx,bb,map.fq,geo.x,geo.b); % fq [Vs]
    end
end

geo.Ar=min(Ar);

RQ=RQ(:,~errorFlag); % filter the unfeasible machines
xRaw=xRaw(~errorFlag);
bRaw=bRaw(~errorFlag);

save('dataSet','dataSet');

if ~isempty(RQ)
    if FEAfixN~=1000 && isempty(isprop(gcp('nocreate'),'NumWorkers'))
        for mot=1:length(xRaw)
            disp(['FEA simulation ' int2str(mot) ' of ' int2str(length(xRaw))])
            [~,~,~,out,~] = FEMMfitness(RQ(:,mot),geo,per,mat,eval_type);
            fdFEA(mot)=out.fd;
            fqFEA(mot)=out.fq;
        end
    else
        parfor mot=1:length(xRaw)
            disp(['FEA simulation ' int2str(mot) ' of ' int2str(length(xRaw))])
            [~,~,~,out,~] = FEMMfitness(RQ(:,mot),geo,per,mat,eval_type);
            fdFEA(mot)=out.fd;
            fqFEA(mot)=out.fq;
        end
    end
else
    disp('FEAfix motor on the x-b plane are not drawable')
    disp('Please correct the x and b ranges to improve the FEA calibration')
    warning('Calibration error')
end


% xRaw=[xRaw(1:2);xRaw(3:4)];
% bRaw=[bRaw(1:2);bRaw(3:4)];
% kdRaw=[kdRaw(1:2);kdRaw(3:4)];
% kqRaw=[kqRaw(1:2);kqRaw(3:4)];
% 
% kd=interp2(xRaw,bRaw,kdRaw,xx,bb);
% kq=interp2(xRaw,bRaw,kqRaw,xx,bb);

kdRaw=fdFEA./fdMod;
kqRaw=fqFEA./fqMod;

kdRaw=kdRaw(~errorFlag);
kqRaw=kqRaw(~errorFlag);

if FEAfixN==1
    kd=kdRaw*ones(size(xx));
    kq=kqRaw*ones(size(xx));
else
    kd=scatteredInterpolant(xRaw',bRaw',kdRaw','linear');
    kq=scatteredInterpolant(xRaw',bRaw',kqRaw','linear');
    
    kd=kd(xx,bb);
    kq=kq(xx,bb);
end

disp('End of FEAfix calibration')

delete('dataSet.mat')

FEAfixOut.kd=kd;
FEAfixOut.kq=kq;
FEAfixOut.xRaw=xRaw;
FEAfixOut.bRaw=bRaw;

if debug
    figure()
    figSetting(12,12)
    xlabel('$x$')
    ylabel('$b$')
    zlabel('$k_d$')
    view(45,30)
    set(gca,'XLim',[min(min(xx)) max(max(xx))],'YLim',[min(min(bb)) max(max(bb))]);
    surf(xx,bb,kd);
    plot3(xRaw,bRaw,kdRaw,'ko')
    
    figure()
    figSetting(12,12)
    xlabel('$x$')
    ylabel('$b$')
    zlabel('$k_d$')
    view(45,30)
    set(gca,'XLim',[min(min(xx)) max(max(xx))],'YLim',[min(min(bb)) max(max(bb))]);
    surf(xx,bb,kq);
    plot3(xRaw,bRaw,kqRaw,'ko')
    
    keyboard
end

timeEnd=toc();
clc
disp(['FEAfix procedure'])
disp(['- number of FEA machines : ' int2str(length(xRaw))])
disp(['- elapsed time           : ' num2str(timeEnd,2) ' s'])

%keyboard

