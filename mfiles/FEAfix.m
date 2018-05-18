function [FEAfixOut] = FEAfix(dataSet,geo,map,FEAfixN)
% 
% [kd,kq] = FEAfix(dataSet,geo,map,setup)
% 

disp('FEAfix calibration...')

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
        xRaw=xx;
        bRaw=bb;
    otherwise
        error('Put a correct number!!!')      
end
kdRaw=ones(size(xRaw));
kqRaw=ones(size(xRaw));
errorFlag=0;


for mot=1:FEAfixN
    disp(['FEA simulation ' int2str(mot) ' of ' int2str(FEAfixN)])
    geo.x=xRaw(mot);
    geo.b=bRaw(mot);

    geo.r = geo.x*geo.R;                                        % rotor radius [mm]
    geo.wt = interp2(xx,bb,map.wt,geo.x,geo.b);                       % tooth width
    geo.wt = round(geo.wt*100)/100;
    geo.lt=interp2(xx,bb,map.lt,geo.x,geo.b);                         % slot length
    geo.lt=round(geo.lt*100)/100;
    geo.x0=geo.R * geo.x /cos(pi/2/geo.p);
    geo.Ar=interp2(xx,bb,map.Ar,geo.x,geo.b);                          % shaft radius [mm]
    geo.Ar=round(geo.Ar*100)/100;
    geo.la=interp2(xx,bb,map.la,geo.x,geo.b);                         % total insulation
    
    if isnan(geo.lt)
        errorFlag=1;
        kdRaw(mot)=1;
        kqRaw(mot)=1;
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

        % current phase angle
        temp_id = interp2(xx,bb,map.id,geo.x,geo.b);         % id [A]
        temp_iq = interp2(xx,bb,map.iq,geo.x,geo.b);         % iq [A]
        temp_fd = interp2(xx,bb,map.Ld,geo.x,geo.b)*temp_id; % fd [Vs]
        temp_fq = interp2(xx,bb,map.Lq,geo.x,geo.b)*temp_iq; % fq [Vs]
        dataSet.GammaPP=round(atand(temp_iq/temp_id)*100)/100;

        % adjourn dataSet
        dataSet.AirGapRadius=round(geo.r*100)/100;
        dataSet.ShaftRadius = round(geo.Ar*100)/100;
        dataSet.ToothLength=geo.lt;
        dataSet.ToothWidth=geo.wt;
        dataSet.ThicknessOfPM = geo.lm;
        dataSet.HCpu=round(geo.hc_pu*100)/100;
        dataSet.DepthOfBarrier=round(geo.dx*100)/100;

        dataSet.currentfilename='mot_FEAfix.mat';
        dataSet.currentpathname='tmp\';

        dataSet.Mesh=2;
        dataSet.Mesh_MOOA=10;

        DrawMachineScript(dataSet,dataSet.currentpathname,dataSet.currentfilename);
        
        eval_type='singt';
        filemot=[pwd '\' dataSet.currentpathname dataSet.currentfilename(1:end-4) '.fem'];
        geo0=geo;
        dataSet0=dataSet;
        per0=per;
        load([filemot(1:end-4) '.mat'])
        per.gamma=dataSet0.GammaPP;
        geo.nsim_singt=4;
        geo.delta_sim_singt=60;
        per.overload=dataSet.CurrLoPP;
        [~,~,~,out,~]=FEMMfitness([],geo,per,mat,eval_type,filemot);
        
        kdRaw(mot)=out.fd/temp_fd;
        kqRaw(mot)=out.fq/temp_fq;
        
        geo=geo0;
        dataSet=dataSet0;
        per=per0;
        
    end
end

if errorFlag
    disp('Some motor on the x-b plane are not drawable')
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

%keyboard

