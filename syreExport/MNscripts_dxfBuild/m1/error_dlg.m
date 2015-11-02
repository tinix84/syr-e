switch error_code(1)
    case 1
        errordlg(['ABORT - il layer ' num2str(error_code(2)) ' sconfina nel successivo'], 'ERROR');
        break
    case 2
        errordlg(['ABORT - il layer 1 è troppo piccolo: aumentare \alpha1 o ridurre hc_1'], 'ERROR');
        break
    case 3
        errordlg(['ABORT - il magnete ' num2str(error_code(2)) ' è troppo piccolo e cozza nel ponticello '], 'ERROR');
        break
    case 4
        errordlg(['ABORT - il layer esterno mi tocca'], 'ERROR');
        break        
end
