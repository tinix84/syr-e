function FemmProblem=draw_lines_archesX(FemmProblem,Mat,gruppo,res)
% Importa_Mat(Mat,gruppo)
[nrig,ncol] = size(Mat);
stop=130;
for ii=1:nrig
%     if (nrig<160)&&(ii==stop)
%         disp('X')
%         writefemmfile(['testZeroX' num2str(ii) '.fem'], FemmProblem)
%         system(['move testZeroX' num2str(ii) '.fem' ' ..'])
%         cd('..')
%         error('OK X')
%     end
    if Mat(ii,ncol)==0
        outernodes = [ Mat(ii,3), Mat(ii,4); Mat(ii,1), Mat(ii,2) ];
        if(norm(outernodes(1,:)-outernodes(2,:))>eps)
            [FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, outernodes(:,1), outernodes(:,2),'InGroup',gruppo);
            [FemmProblem] = addsegments_mfemm(FemmProblem, nodeids(1), nodeids(2),...
                'MaxSideLength',res,'Hidden',0,'InGroup',gruppo,'BoundaryMarker','None');
        end
        %mi_drawline(Mat(ii,3),Mat(ii,4),Mat(ii,1),Mat(ii,2));
        %mi_selectsegment(mean([Mat(ii,1) Mat(ii,3)]),mean([Mat(ii,2) Mat(ii,4)]));
        % 21-06-2010 BB
        % La dimensione degli elementi della mesh non e' piu' impostata in
        % automatico dal programma, ma fissata in base al valore di input
        % res
        % mi_setsegmentprop('None', 0, 1, 0, gruppo);
        %         mi_setsegmentprop('None', res, 0, 0, gruppo);
        %         mi_clearselected
        %         mi_selectnode(Mat(ii,1),Mat(ii,2));
        %         mi_selectnode(Mat(ii,3),Mat(ii,4));
        %         mi_setnodeprop('None',gruppo);
        %         mi_clearselected
        %         keyboard
    else
        [outernodes,maxsegdeg,~,~,ang]=Disegna_ArcoX(Mat(ii,:),res);
        [FemmProblem, ~, nodeids] = addnodes_mfemm(FemmProblem, ...
            outernodes(:,1), outernodes(:,2),'InGroup',gruppo);
        [FemmProblem, rcsegind] = addarcsegments_mfemm(FemmProblem, nodeids(1), ...
            nodeids(2), ang, 'MaxSegDegrees', maxsegdeg, ...
            'InGroup', gruppo, 'BoundaryMarker','None');
        
        %         [x_temp,y_temp]=pol2cart((ang1+0.5*ang)*pi/180,raggio);
        %         x=x_temp+Mat(ii,1);
        %         y=y_temp+Mat(ii,2);
        %         % mi_selectarcsegment( Mat(ii,5)*0.01+Mat(ii,3)*0.99 , Mat(ii,6)*0.01+Mat(ii,4)*0.99 );
        %         % mi_setarcsegmentprop(res, 'None', 0, gruppo);
        %         mi_selectarcsegment(x,y);
        %         mi_setarcsegmentprop(maxsegdeg, 'None', 0, gruppo);
        %         mi_clearselected
        %         mi_selectnode(Mat(ii,5),Mat(ii,6));
        %         mi_selectnode(Mat(ii,3),Mat(ii,4));
        %         mi_setnodeprop('None',gruppo);
        %         mi_clearselected
        %        keyboard
    end
end
% writefemmfile('testZeroX.fem', FemmProblem)
% system(['move testZeroX.fem ..'])
% cd('..')
% error('OK X')





