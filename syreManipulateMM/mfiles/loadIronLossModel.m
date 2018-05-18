function [IronLossModel]=loadIronLossModel(filename)

if strcmp(filename,'0')
    IronLossModel.type='0';
else
    load(filename)
    if exist('Pfes_c','var') % model from MagNet maps in idiq
        IronLossModel.type='map';
        %IronLossModel.IdMax=max(max(Id));
        %IronLossModel.IdMin=min(min(Id));
        %IronLossModel.IqMax=max(max(Iq));
        %IronLossModel.IqMin=min(min(Iq));
        IronLossModel.Pfe_h=Pfes_h+Pfer_h;
        IronLossModel.Pfe_c=Pfes_c+Pfer_c;
        IronLossModel.Ppm=Ppm;
        IronLossModel.n0=velDim;
        IronLossModel.Id=Id;
        IronLossModel.Iq=Iq;
        
        prompt={'Pole pair number',...
                'Exponent for eddy current loss', ...
                'Exponent for hysteresis loss', ...
                'Exponent for PM loss', ...
                'Segmentation factor for PM loss'};
        name='Input';
        numlines=1;
        defaultanswer={'2','2','1.2','2','1'};
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        
        p=eval(setup{1});
        IronLossModel.f0=IronLossModel.n0/60*p;
        IronLossModel.expC=eval(setup{2});
        IronLossModel.expH=eval(setup{3});
        IronLossModel.expPM=eval(setup{4});
        IronLossModel.segPM=eval(setup{5});
        
    elseif exist('IronLossModel','var') % model from fit
        if ~isfield(IronLossModel,'type')
            if isfield(IronLossModel,'param')
                IronLossModel.type='fitIM1';
            else
                error('Model not valid')
            end
        end
    elseif exist('out','var') %single point, to verify
        IronLossModel.type='point';
        IronLossModel.Pfe0=out.Pfes_h+out.Pfes_c+out.Pfer_h+out.Pfer_c;
        IronLossModel.F0=abs(out.Fd+j*out.Fq);
        IronLossModel.n0=out.veldim;
        
        prompt={'Pole pair number',...
                'Exponent for flux linkage', ...
                'Exponent for frequency'};
        name='Input';
        numlines=1;
        defaultanswer={'2','2','1.8'};
        setup=inputdlg(prompt,name,numlines,defaultanswer);
        
        p=eval(setup{1});
        IronLossModel.f0=IronLossModel.n0*p/60;
        IronLossModel.expFlux=eval(setup{2});
        IronLossModel.expFreq=eval(setup{3});
    end
end
        