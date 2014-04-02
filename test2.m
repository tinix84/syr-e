if matlabpool('size') == 0
 %matlabpool local
 matlabpool
end

try 
    parfor i = 1:4
        femmHandle = actxserver('femm.ActiveFEMM');
        disp('Server created');
        
        invoke(femmHandle, 'mlab2femm', ['create(' , num(0) , ')' ]); % create new document
        invoke(femmHandle, 'mlab2femm', ['mi_addnode(' , numc(0), num(0), ')']); % invoke some dummy functions
        invoke(femmHandle, 'mlab2femm', ['mi_addnode(' , numc(10), num(10), ')']);
        
        for j = 1:100
            invoke(femmHandle, 'mlab2femm', ['mi_addsegment(' , numc(0) , numc(0) , numc(10) , num(10) , ')']); 
            invoke(femmHandle, 'mlab2femm', ['mi_selectsegment(' , numc(0), num(0), ')']); 
            invoke(femmHandle, 'mlab2femm', ['mi_deleteselected()']); 
        end
        
        delete(femmHandle); 
        disp('Server released');
    end
    disp('Finished');
catch exception
    disp('Fatal Error'); 
    disp(exception.message);
    disp(exception.identifier);
end