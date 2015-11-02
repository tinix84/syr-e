%% tegoli_in_SmCo_ 2mag x layer

Smag=[5.7,7.5,9.9]; %lunghezza magneti

        if (jj==1)
            mi_drawline(X32(jj),Smag(jj),X42(jj),Smag(jj));
            mi_drawline(X32(jj),-Smag(jj),X42(jj),-Smag(jj));
            xa1=X32(jj)-(X32(jj)-X42(jj))/2;
            ya1=Y32(jj)+hc(jj)/4;
            mi_addblocklabel(xa1,ya1);mi_selectlabel(xa1,ya1);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;
            mi_addblocklabel(xa1,-ya1); mi_selectlabel(xa1,-ya1);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;
            
        elseif (jj==2)

            mi_drawline(x32,Smag(jj),x42,Smag(jj));
            mi_drawline(x32,-Smag(jj),x42,-Smag(jj));
            
            nax_est=X3(jj);
            nay_est=Y3(jj)+hc(jj)/2*cos(pi/4);
            mi_addblocklabel(nax_est,nay_est);
            mi_selectlabel(nax_est,nay_est);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;
            mi_addblocklabel(nax_est,-nay_est);
            mi_selectlabel(nax_est,-nay_est);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;

            
        else
            mi_drawline(x32,Smag(jj),x42,Smag(jj));
            mi_drawline(x32,-Smag(jj),x42,-Smag(jj));
            
            % aria vicino ai ponticelli nei layer radiali
            nax_est=X3(jj);
            nay_est=Y3(jj)+hc(jj)/2*cos(pi/4);
            mi_addblocklabel(nax_est,nay_est);
             mi_selectlabel(nax_est,nay_est);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;
            mi_addblocklabel(nax_est,-nay_est);
            mi_selectlabel(nax_est,-nay_est);
            mi_setblockprop('Air', 1, 0, 'None', 0, group0, 1);
            mi_clearselected;
           
        end
