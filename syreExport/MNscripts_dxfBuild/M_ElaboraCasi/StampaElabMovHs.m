
mkdir([FilePath 'fig']);

print(1,'-dpdf','-r300',[FilePath 'fig\01_Coppia-a'])
print(2,'-dpdf','-r300',[FilePath 'fig\02_Flussi-a'])
print(3,'-dpdf','-r300',[FilePath 'fig\03_EMF-a'])
print(4,'-dpdf','-r300',[FilePath 'fig\04_FluxDQ-a'])
print(5,'-dpdf','-r300',[FilePath 'fig\05_Coppia-FFT'])
print(6,'-dpdf','-r300',[FilePath 'fig\06_EMF-FFT'])

saveas(1,[FilePath 'fig\01_Coppia-a'])
saveas(2,[FilePath 'fig\02_Flussi-a'])
saveas(3,[FilePath 'fig\03_EMF-a'])
saveas(4,[FilePath 'fig\04_FluxDQ-a'])
saveas(5,[FilePath 'fig\05_Coppia-FFT'])
saveas(6,[FilePath 'fig\06_EMF-FFT'])


% % Costruisce un file unico
% switch computer
%  case {'GLNXA64','GLNX86'}
%      StringaComando=['!gs -q -dNOPAUSE -dSAFER -dBATCH -sOutputFile=' FilePath Mac.MachineName '.pdf -sDEVICE=pdfwrite ' FilePath '*.pdf'];
%      eval(StringaComando)
%      StringaComando=['!rm ' FilePath '0*.pdf'];
%      eval(StringaComando) 
%  case 'PCWIN'
%      % pdftk: merge test1.pdf and test2.pdf into output.pdf
%      % !pdftk test1.pdf test2.pdf cat output output.pdf
%      ElencoNomiFile = ['01_Coppia-a.pdf 02_Flussi-a.pdf 03_EMF-a.pdf 04_FluxDQ-a.pdf 05_Coppia-t.pdf 06_Coppia-s.pdf 07_EMF-t.pdf 08_EMF-s.pdf']
%      ElencoNomiFile = [FilePath '*.pdf'];
%      StringaComando=['!pdftk ' ElencoNomiFile ' cat output ' FilePath 'output.pdf'];
%      eval(StringaComando);
%      clear ElencoNomiFile StringaComando
% end