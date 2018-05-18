function Pfe=calcIronLoss(IronLossModel,fdfq,FreqElet)

Id=fdfq.Id;
Iq=fdfq.Iq;
Fd=fdfq.Fd;
Fq=fdfq.Fq;

if strcmp(IronLossModel.type,'0')
    Pfe=zeros(size(Id));
elseif strcmp(IronLossModel.type,'map')
    IdMap=IronLossModel.Id;
    IqMap=IronLossModel.Iq;
    
    %if max(max(Id))>max(max(IdMap))
    %    flag=1;
    %elseif min(min(Id))<min(min(IdMap))
    %    flag=1;
    %elseif max(max(Iq))>max(max(IqMap))
    %    flag=1;
    %elseif min(min(Iq))<min(min(IqMap))
    %    flag=1;
    %else
    %    flag=0;
    %end
    
    Pfeh = interp2(IdMap,IqMap,IronLossModel.Pfe_h,Id,Iq,'cubic',1e50);
    Pfec = interp2(IdMap,IqMap,IronLossModel.Pfe_c,Id,Iq,'cubic',1e50);
    Ppm  = interp2(IdMap,IqMap,IronLossModel.Ppm,Id,Iq,'cubic',1e50);

    Pfeh = Pfeh*(FreqElet./IronLossModel.f0).^IronLossModel.expH;
    Pfec = Pfec*(FreqElet./IronLossModel.f0).^IronLossModel.expC;
    Ppm  = IronLossModel.segPM*Ppm*(FreqElet./IronLossModel.f0).^IronLossModel.expPM;

    Pfe=Pfeh+Pfec+Ppm;
    
    %if flag
    %    Pfe=1e50;
    %else
    %   Pfeh = interp2(IdMap,IqMap,IronLossModel.Pfe_h,Id,Iq);
    %   Pfec = interp2(IdMap,IqMap,IronLossModel.Pfe_c,Id,Iq);
    %    Ppm  = interp2(IdMap,IqMap,IronLossModel.Ppm,Id,Iq);
    %   
    %   Pfeh = Pfeh*(FreqElet./IronLossModel.f0).^IronLossModel.expH;
    %    Pfec = Pfec*(FreqElet./IronLossModel.f0).^IronLossModel.expC;
    %    Ppm  = Ppm*(FreqElet./IronLossModel.f0).^IronLossModel.expPM;
    %    
    %    Pfe=Pfeh+Pfec+Ppm;
    %end
elseif strcmp(IronLossModel.type,'fitIM1')
    a  = IronLossModel.param.a;
    b  = IronLossModel.param.b;
    k1 = IronLossModel.param.k1;
    k2 = IronLossModel.param.k2;
    k3 = IronLossModel.param.k3;
    f0 = IronLossModel.param.f0;
    c1 = IronLossModel.param.c1;
    c2 = IronLossModel.param.c2;
    c3 = IronLossModel.param.c3;

    Pfe=k1*(fs/f0).^c1.*Fs.^a+k2*(fs/f0).^c2.*Iq.^b-k3*(fs/f0).^c3.*Fs.*Iq;
elseif strcmp(IronLossModel.type,'point')
    Pfe0=IronLossModel.Pfe0;
    F0=IronLossModel.F0;
    f0=IronLossModel.f0;
    expFlux=IronLossModel.expFlux;
    expFreq=IronLossModel.expFreq;
    
    Fabs=abs(Fd+j*Fq);
    
    Pfe=Pfe0*(Fabs./F0).^expFlux*(FreqElet./f0).^expFreq;
elseif strcmp(IronLossModel.type,'fitExpSyRMTPA')
    k1=IronLossModel.param.k1;
    k2=IronLossModel.param.k2;
    k3=IronLossModel.param.k3;
    e1=IronLossModel.param.e1;
    e2=IronLossModel.param.e2;
    
    Iabs=abs(Id+j*Iq);
    Pfe=k1*FreqElet.^e1+k2*Iabs.^e2+k3*FreqElet.*Iabs;
    Pfe(Pfe<0)=0;
end




