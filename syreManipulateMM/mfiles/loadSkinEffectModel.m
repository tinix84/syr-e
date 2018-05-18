function SkinEffModel=loadSkinEffectModel(filename)

if strcmp(filename,'0')
    SkinEffModel.type='0';
else
    load(filename)
    f=results.f(1,:);
    k=results.k(1,:);
    [p,s] = polyfit(f,k,7);
    
    SkinEffModel.type='interp';
    SkinEffModel.f=f;
    SkinEffModel.k=k;
    SkinEffModel.p=p;
    SkinEffModel.s=s;
    SkinEffModel.n=7;
end
        