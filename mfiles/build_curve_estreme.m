%% 26 12 2011 - GP
if exist('Ld')
    %% ldlq_idiq_n128 - no need for 256x256 INTERPOLATION
    if axes_type == 'PM'
        
        FluxDmin = Ld(1,:);
        
        %     [value, index] = min(abs(Iq(:,1)))
        %     FluxDo = Ld(index,:);
        FluxDo = Ld(1,:);
        FluxDmax = Ld(end,:);
        [value, index] = min(abs(Id(1,:)))
        FluxQo = Lq(:,index);
        FluxQmax = Lq(:,1);
    else
        FluxDmin = Ld(1,:);
        
        [value, index] = min(abs(Iq(:,1)))
        FluxDo = Ld(index,:);
        FluxDmax = Ld(end,:);
        FluxQo = Lq(:,1);
        FluxQmax = Lq(:,end);
    end
    Id = Id(1,:);
    Iq = Iq(:,1);
    
    n_interp = 0;   % number of extra interpolations
    
elseif or(exist('F_map'),exist('out'))
    %% MAGNET and MOGA style (F_map)
    if not(exist('F_map'))
        F_map = out;
    end
    if (ndims(F_map.Fd) == 2)
        %% no skew
        if size(F_map.Fd,1) == 4
            % curve estreme
            FluxDmin = F_map.Fd(1,:);
            FluxDo = F_map.Fd(1,:);
            FluxDmax = F_map.Fd(2,:);
            FluxQo = F_map.Fq(3,:);
            FluxQmax = F_map.Fq(4,:);
            Id = F_map.Id(1,:);
            Iq = F_map.Iq(3,:);
            
            n_interp = 2;   % number of extra interpolations
            
        else
            
            % matrice di punti
            FluxDmin = F_map.Fd(1,:);
            FluxDo = F_map.Fd(1,:);
            FluxDmax = F_map.Fd(end,:);
            FluxQo = F_map.Fq(:,1);
            FluxQmax = F_map.Fq(:,end);
            Id = F_map.Id(1,:);
            Iq = F_map.Iq(:,1);
            
            n_interp = 1;   % number of extra interpolations
            
        end
    elseif (ndims(F_map.Fd) == 3)
        %% skew
        FluxDmin = [];
        FluxDo = mean(F_map.Fd(1,:,:),3);
        FluxDmax = mean(F_map.Fd(2,:,:),3);
        FluxQo = mean(F_map.Fq(3,:,:),3);
        FluxQmax = mean(F_map.Fq(4,:,:),3);
        Id = F_map.Id(1,:,1);
        Iq = F_map.Iq(3,:,1);
        
        n_interp = 2;   % number of extra interpolations
        
    else
        warndlg('wait !!')
    end
elseif exist('i2flux')
    %% pastorelli EXP
    FluxDmin = i2flux.fluxd(1,:);
    ind = find(i2flux.iq_mLd>=0,1,'first');
    FluxDo = i2flux.fluxd(ind,:);
    FluxDmax = i2flux.fluxd(end,:);
    FluxQo = i2flux.fluxq(:,1);
    FluxQmax = i2flux.fluxq(:,end);
    Id = i2flux.id;
    Iq = i2flux.iq;
    
    [F_map.Id,F_map.Iq] = meshgrid(Id,Iq);
    F_map.Fd = i2flux.fluxd; F_map.Fq = i2flux.fluxq;
    
    n_interp = 0;   % number of extra interpolations
    
else
    %% other formats
    if exist('FDtot')   % MappaSally ...
        FluxDmin = FDtot(:,1);
        FluxDo = FDtot(:,IQtot == 0);
        FluxDmax = FDtot(:,end);
        FluxQo = FQtot(1,:);
        FluxQmax = FQtot(end,:);
        Id = IDtot;
        Iq = IQtot;
        
        n_interp = 0;   % number of extra interpolations
        
    elseif exist('F')   %~isempty('F')      % mappa magnet CE1SM
        Id =        F.values(1:10,1)';
        Iq =        F.values(21:30,2)';
        FluxDmin =  [];
        FluxDo  =   F.values(1:10,3)';
        FluxDmax =  F.values(11:20,3)';
        FluxQoM  =   F.values(21:30,4)';
        FluxQmax = 	F.values(31:40,4)';
        
        n_interp = 0;   % number of extra interpolations
        
    elseif exist('Ld_Exp')   %~isempty('F')      % mappa magnet CE1SM
        
        %% inverter E3
        FluxDmin = Ld_Exp(1,:);
        ind = find(Iq_Exp>=0,1,'first');
        FluxDo = Ld_Exp(ind,:);
        FluxDmax = Ld_Exp(end,:);
        FluxQo = Lq_Exp(:,1);
        FluxQmax = Lq_Exp(:,end);
        Id = Id_Exp(1,:);
        Iq = Iq_Exp(:,1);
        
        [F_map.Id,F_map.Iq] = meshgrid(Id,Iq);
        F_map.Fd = Ld_Exp; F_map.Fq = Lq_Exp;
        
        n_interp = 1;   % number of extra interpolations
        
    elseif exist('Output')   %~isempty('F')      % mappa magnet CE1SM
        
        FluxDmin = [];
        FluxDo = [];
        FluxDmax = [];
        FluxQo = [];
        FluxQmax = [];
        
        Id = Output.id(1,:);
        Iq = Output.iq(:,1);
        
        F_map.Id = Output.id;
        F_map.Iq = Output.iq;
        if isfield(Output.MagModelId,'LambdaVRef_bk')
            F_map.Fd = 0.5 * ((Output.MagModelId.LambdaVRef.d) + (Output.MagModelId.LambdaVRef_bk.d));
            F_map.Fq = 0.5 * (abs(Output.MagModelId.LambdaVRef.q) + abs(Output.MagModelId.LambdaVRef_bk.q));
            F_map.T = 0.5 * (abs(Output.MagModelId.Torque) + abs(Output.MagModelId.Torque_bk));
        else
            F_map.Fd = Output.MagModelId.LambdaVRef.d;
            F_map.Fq = Output.MagModelId.LambdaVRef.q;
        end
        
        n_interp = 1;   % number of extra interpolations
        
        
    end
end

if ~isempty(FluxDo)
    %% end-winding inductances
    FluxDmin =  FluxDmin + Lld * Id;
    FluxDo  =   FluxDo + Lld * Id;
    FluxDmax =  FluxDmax + Lld * Id;
    FluxQo =  FluxQo + Llq * Iq;
    FluxQmax = 	FluxQmax + Llq * Iq;
    
    figure, hold on
    if ~isempty(FluxDmin)
        plot(Id,FluxDmin,'-','Color',[0 0.5 0])
    end
    plot(Id,FluxDo,'b-','LineWidth',[2]),grid on,
    plot(Id,FluxDmax,'-','Color',[0 0.5 0],'LineWidth',[2])
    plot(Iq,FluxQo,'b-','LineWidth',[2])
    plot(Iq,FluxQmax,'-','Color',[0 0.5 0],'LineWidth',[2])
    xlabel('Apk'), ylabel('Vs')
    
    
    % descr = input('motor name ??')
    % comment = input('comment: ')
    % title([descr ' - magnetic model - ' comment]);
    title([motor_name ' - magnetic model']);
    saveas(gcf,[pathname1 'CurveEstreme_' motor_name '.fig']);
    save([pathname1 'CurveEstreme.mat'],'Id','Iq','FluxDo','FluxDmax','FluxQo','FluxQmax');
end
