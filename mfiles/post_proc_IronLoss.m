% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

% load phase flux linkages
temp_out = mo_getcircuitproperties('fase1');
temp_out = temp_out - mo_getcircuitproperties('fase1n');
f1 = temp_out(3) * 2 * p/ps;
temp_out = mo_getcircuitproperties('fase2');
temp_out = temp_out - mo_getcircuitproperties('fase2n');
f2 = temp_out(3) * 2 * p/ps;
temp_out = mo_getcircuitproperties('fase3');
temp_out = temp_out - mo_getcircuitproperties('fase3n');
f3 = temp_out(3) * 2 * p/ps;

% evaluate torque
% T1 - from the innermost integration line (rot + gap/6)
x = r + gap*1/6;
ang0 = th_m; ang1 = gradi_da_sim + th_m;
[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);

T1 = mo_lineintegral(4);
T1 = T1(1) * 2 * p/ps;
mo_clearcontour();

% T2 - from the outermost integration line (stat - gap/6)
x = r + gap*5/6;
ang0 = -pc; ang1 = gradi_da_sim-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);
T2 = mo_lineintegral(4);
T2 = T2(1) * 2 * p/ps;
mo_clearcontour();

% T3 - from an intermediate line (rot + gap/2)
x = r + gap*1/2;
ang0 = -pc; ang1 = gradi_da_sim-pc;

[x1,y1] = rot_point(x,0,ang0*pi/180);
[x2,y2] = rot_point(x,0,ang1*pi/180);
mo_addcontour(x1,y1);
mo_addcontour(x2,y2);
mo_bendcontour(gradi_da_sim,0.5);
T3 = mo_lineintegral(4);
T3 = T3(1) * 2 * p/ps;
mo_clearcontour();

% dq flux linkaged
fdq = abc2dq(f1,f2,f3,th(jj)*pi/180);

% string of SOL

sol = [th(jj) id iq fdq(1) fdq(2) mean([T1,T2,T3])];

% % Store first simulation status, including sampling coordinate, group
% number, saved in Matrix
if jj == 1
    EleNo = mo_numelements;               % Number of mesh elements
    pos = zeros(EleNo,1);                 % Matrix that will hold the mesh elements centroid coordinates as complex number
    area = zeros(EleNo,1);                 % Matrix that will hold the mesh elements area
    groNo = zeros(EleNo,1);               % Matrix that will hold the mesh elements group number
    for i = 1:EleNo
        elm = mo_getelement(i);
        pos(i) = elm(4)+j*elm(5);
        area(i) = elm(6);
        groNo(i) = elm(7);
    end
  
    %% find stator iron elements
    [line_Sta, row] = find(groNo==12);      % Stator Iron
    pos_Sta = pos(line_Sta);
    area_Sta = area(line_Sta);
    gro_Sta = groNo(line_Sta);
    %% select 1/3 stator iron elements
    Sta = find(angle(pos_Sta) < ((Qs/3*2-1)*pc*pi/180));
    pos_Sta = pos_Sta(Sta);
    area_Sta = area_Sta(Sta);
    gro_Sta = gro_Sta(Sta);
    
    %% find rotor iron element
    [line_Rot, row] = find(groNo==22);
    pos_Rot = pos(line_Rot);
    area_Rot = area(line_Rot);
    gro_Rot = groNo(line_Rot);
    
    %% Matrix to store data
    areaIron = [area_Sta; area_Rot];
    posIron = [pos_Sta; pos_Rot];
    groNoIron = [gro_Sta; gro_Rot];
    %% delete redundant data, just keep iron element data
    IronNo = size(groNoIron,1);
    fIron =zeros(IronNo,nsim);            % Matrix that save flux density of iron element
end

RotPos = exp(j*th_m*pi/180);              % It is a parameter that considers the rotor position

for i = 1:IronNo
    if (groNoIron(i) == 22)                       % rotor iron elements
        Pos_Rot = posIron(i)*RotPos;              % get element position
        %% axis transform
        fIron(i,jj) = (mo_getb(real(Pos_Rot),imag(Pos_Rot))*[1;j]);
        Bd = real(fIron(i,jj)) * cosd(th_m) + imag(fIron(i,jj)) * sind(th_m);
        Bq = imag(fIron(i,jj)) * cosd(th_m) - real(fIron(i,jj)) * sind(th_m);
        fIron(i,jj) = [Bd,Bq] * [1;j];
    elseif (groNoIron(i) == 12)                   % stator iron elements
        if ceil(th_m) < 60/p    % for stator, just 60 ele degree rotation data is needed
%%          first one third of stator
            Pos_Sta = posIron(i);
            fIron(i,jj) = (mo_getb(real(Pos_Sta),imag(Pos_Sta))*[1;j]);
%%          second one third stator iron area
            Pos_Sta1 = Pos_Sta*exp(j*(Qs/3*2*pc)*pi/180);
            fIron(i,jj+nsim) = (mo_getb(real(Pos_Sta1),imag(Pos_Sta1))*[1;j]);
%%          last one third stator iron area
            Pos_Sta2 = Pos_Sta1*exp(j*(Qs/3*2*pc)*pi/180);
            fIron(i,jj+nsim*2) = (mo_getb(real(Pos_Sta2),imag(Pos_Sta2))*[1;j]);
        end
    end
end

if jj == nsim
    fluxdens = [posIron areaIron groNoIron fIron];      % Saving all the info necessary to calculate Iron Losses
end