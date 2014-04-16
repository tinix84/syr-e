% draw_motor_in_FEMM.m - 24 08 09 - GMP
% builds the .fem motor model with the rotor in position zero

% input: geo and mat
% output: mot0.fem (th_m = 0, i123 = [0,0,0])

ns = geo.ns;
p = geo.p;

fem = dimMesh(geo,eval_type);
% fem.res_traf=0.5;
% fem.res=6;
filename = 'mot0.fem';
%opendocument('empty_case.fem');
%mi_probdef(0,'millimeters','planar',1e-8,geo.l,25);

% calc winding factor (kavv) and rotor offset respect to position zero
% (phase1_offset)
[kavv, phase1_offset] = calcola_kavv_th0(geo.avv,geo.ns*geo.p,size(geo.avv,1),geo.p);

% rotor offset angle
th_m0 = 0;                          % [mec deg]
% offset angle: for coordinate transformations
geo.th0 = th_m0*p - phase1_offset;     % [elt deg]
%%
%% evaluation of the number of poles to simulate and corresponding periodicity
%%
Q=geo.ns*geo.p;
t=gcd(ns*geo.p,geo.p);  % periodicity
if ((6*t/Q)>1)
    ps=2*p/t;
    Qs=Q/t;
else
    ps=p/t;
    Qs=Q/2/t;
end
%assign to variable geo, number of poles and slots simulated.
geo.Qs=Qs;
geo.ps=ps;
% Boundary tipology
if (rem(geo.ps,2)==0)
    periodicity=4;
else
    periodicity=5;
end

%% definition of the boundary conditions
% inner and outer circles
[FemmProblem, boundindA0, boundnameA0]=addboundaryprop_mfemm(FemmProblem,'A=0',0);
%mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);
% Anit-Periodicity (2 x rotor + 2 x stator + 3 x airgap + 1 x AP move, which is the sliding contour)
[FemmProblem, boundindAPr1, boundnameAPr1]=addboundaryprop_mfemm(FemmProblem,'APr1',periodicity);
%mi_addboundprop('APr1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPr2, boundnameAPr2]=addboundaryprop_mfemm(FemmProblem,'APr2',periodicity);
%mi_addboundprop('APr2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPg1, boundnameAPg1]=addboundaryprop_mfemm(FemmProblem,'APg1',periodicity);
%mi_addboundprop('APg1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPg2, boundnameAPg2]=addboundaryprop_mfemm(FemmProblem,'APg2',periodicity);
%mi_addboundprop('APg2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPg3, boundnameAPg3]=addboundaryprop_mfemm(FemmProblem,'APg3',periodicity);
%mi_addboundprop('APg3', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPmove, boundnameAPmove]=addboundaryprop_mfemm(FemmProblem,'APmove',periodicity);
%mi_addboundprop('APmove', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPs1, boundnameAPs1]=addboundaryprop_mfemm(FemmProblem,'APs1',periodicity);
%mi_addboundprop('APs1', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);
[FemmProblem, boundindAPs2, boundnameAPs2]=addboundaryprop_mfemm(FemmProblem,'APs2',periodicity);
%mi_addboundprop('APs2', 0, 0, 0, 0, 0, 0, 0, 0, periodicity);

%% ROTOR
% build the matrixes which describe the rotor
ROTmatr;
save CaseSim.mat;   % 2014/02/25 MG in the final version this save MUST be eliminated (necessary becouse if the algorith end for error there you find the last case analysed.
BLKLABELS.rotore=BLKLABELSrot;
% draws lines and arces
FemmProblem=draw_lines_archesX(FemmProblem,rotore2,2,fem.res);
% assigns the block labels
%plotfemmproblem(FemmProblem)
FemmProblem=assign_block_prop_rotX(FemmProblem,BLKLABELS,fem,2);
% boundary conditions
for ii=1:2
    %mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.Segments(id).BoundaryMarker = boundnameAPr1;
        %mi_setsegmentprop('APr1', 0, 1, 0, 2);
        %mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        [id] = findarcsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.ArcSegments(id).BoundaryMarker = boundnameA0;
%         mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
%         mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
%         mi_clearselected;
    end
end
% Condizioni al contorno di rotore ferro lamierino
for ii=3:4
    
    %mi_selectsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2));
    if (BLKLABELSrot.boundary(ii,3)==10)
        [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        FemmProblem.Segments(id).BoundaryMarker = boundnameAPr2;
        %mi_setsegmentprop('APr2', 0, 1, 0, 2);
        %mi_clearselected;
    elseif(BLKLABELSrot.boundary(ii,3)==0)
        [id] = findarcsegment_mfemm(FemmProblem, [BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2)]);
        %mi_selectarcsegment(BLKLABELSrot.boundary(ii,1),BLKLABELSrot.boundary(ii,2))
        FemmProblem.ArcSegments(id).BoundaryMarker = boundnameA0;
        %mi_setarcsegmentprop(fem.res, 'A=0', 0, 2);
        %mi_clearselected;
    end
end

%% STATOR
% builds the matrixed which describe the stator
STATmatr;
BLKLABELS.statore=BLKLABELSstat;
% draws lines and arces
FemmProblem=draw_lines_archesX(FemmProblem,statore,1,fem.res);
% assigns the block labels
FemmProblem=assign_block_prop_statX(FemmProblem,BLKLABELS,geo,fem,1); % assegna materiali;
% boundary conditions
BLKLABELSstat=BLKLABELS.statore;
for ii=1:size(BLKLABELSstat.boundary,1)
    [id, xycoords] = findsegment_mfemm(FemmProblem, [BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2)]);
    %mi_selectsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
    if (BLKLABELSstat.boundary(ii,3)==10)
        FemmProblem.Segments(id).BoundaryMarker = boundnameAPs1;
        %mi_setsegmentprop('APs1', 0, 1, 0, 1);
        %mi_clearselected;
    elseif(BLKLABELSstat.boundary(ii,3)==0)
        [id] = findarcsegment_mfemm(FemmProblem, [BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2)]);
        FemmProblem.ArcSegments(id).BoundaryMarker = boundnameA0;
%         mi_selectarcsegment(BLKLABELSstat.boundary(ii,1),BLKLABELSstat.boundary(ii,2));
%         mi_setarcsegmentprop(fem.res, 'A=0', 0, 1);
%         mi_clearselected;
    end
end
% keyboard
%% airgap (group 20)
    FemmProblem=AirGapBuildX(FemmProblem,Qs,ps,geo.p,geo.g,360/(ns*geo.p)/2,geo.xr,fem.res_traf,1,2,boundnameAPg1,boundnameAPg2,boundnameAPg3);
    FemmProblem=draw_airgap_arc_with_meshX(FemmProblem,geo,th_m0,fem,boundnameAPmove);

geo.fem=fem;
%mi_saveas(filename); % saves the file with name ’filename’.

