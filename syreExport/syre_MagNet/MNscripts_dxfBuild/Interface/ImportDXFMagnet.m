function ImportDXFMagnet(h,dxfName)

% ImportDXFMagnet.m [v1.00.00 (03-05-2013)]
% Imports a DXF
% =========================================================================
% ImportDXFMagnet(h,dxfName)
% Input:
%          - h: current Magnet's data structure
%          - dxfName: path and name of the DXF file to be imported
% =========================================================================
mh = h.magnetHandler;
invoke(mh,'processCommand',['CALL getDocument().getView().importDXF("',dxfName,'")']);