% adapts 
function fem = dimMesh(geo)

global eval_type 

res_max=1/3*geo.g;
% La dimensione degli elementi della mesh al traferro deve essere un
% multiplo (o sottomultiplo) del passo di simulazione.
if (eval_type == 'MO_OA')
    fem.res_traf=1/geo.p;
    fem.res=geo.K_mesh_MOOA*fem.res_traf;
else
    res_multiplo=(geo.delta_sim_singt/(geo.nsim_singt-1)*180/pi)*...
        (geo.xr+0.5*geo.g);
    if res_multiplo<1.5*res_max
        fem.res_traf=res_multiplo;
    else
        % Se res_multiplo e' maggiore del massimo (res_max), impongo che il
        % valore di risoluzione scelto per il traferro sia un sottomultiplo di
        % res_multiplo
        cerco_res=res_multiplo;
        divido=2;
        while cerco_res>1.5*res_max
            cerco_res=res_multiplo/divido;
            divido=divido+1;
        end
        fem.res_traf=cerco_res;
    end
    % - Resto della macchina
    fem.res=geo.K_mesh*fem.res_traf;
end

