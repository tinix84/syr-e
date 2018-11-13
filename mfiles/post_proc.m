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
for ii=0:(n3phase-1) %AS
    temp_out = mo_getcircuitproperties(phase_name{3*ii+1});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+1});
    f(3*ii+1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
    temp_out = mo_getcircuitproperties(phase_name{3*ii+2});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+2});
    f(3*ii+2) = temp_out(3) * 2 * p/ps;
    temp_out = mo_getcircuitproperties(phase_name{3*ii+3});
    temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+3});
    f(3*ii+3) = temp_out(3) * 2 * p/ps;
end

% dq flux linkage
for ik=0:(n3phase-1) %AS
    fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),th(jj)*pi/180,n3phase,ik);
    fd_temp(ik+1,jj)=fdq(1);
    fq_temp(ik+1,jj)=fdq(2);
end

fd=mean(fd_temp(:,jj));
fq=mean(fq_temp(:,jj));

% block evaluation of torque and force
for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
    xB=geo.BLKLABELS.rotore.xy(ii,1);
    yB=geo.BLKLABELS.rotore.xy(ii,2);
    [xB,yB]=rot_point(xB,yB,th_m*pi/180);
    mo_selectblock(xB,yB);
end

Tblock=mo_blockintegral(22)*2*p/ps;
mo_clearblock;

% Calcolo del Volume dei magneti VolPM - rev.Gallo 14/03/2018
if jj==1 %viene calcolato solo alla prima simulazione (rotore in posizione di partenza)
    flagPM=0;
    for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
        if geo.BLKLABELS.rotore.xy(ii,3)== 6 %lettura codice del materiale per individuare regioni di PM presenti nella struttura di rotore
            xA=geo.BLKLABELS.rotore.xy(ii,1);
            yA=geo.BLKLABELS.rotore.xy(ii,2);
            [xA,yA]=rot_point(xA,yA,th_m*pi/180); %rotazione delle coordinate di random position offset
            mo_selectblock(xA,yA);
            flagPM=1;
        end
    end
    if flagPM %flag per capire se area del magnete rettangolare è presente o no
        VolPM=(2*geo.p*mo_blockintegral(10))/geo.ps; %Calcolo Volume magnete totale nel rotore [m3]
    else
        VolPM=0;
    end
    
    mo_clearblock();
end

if(0)
    % % Store first simulation status, including sampling coordinate, group
    % number, saved in Matrix
    if jj == 1
        EleNo = mo_numelements;               % Number of mesh elements
        pos = zeros(EleNo,1);                 % Matrix that will hold the mesh elements centroid coordinates as complex number
        area = zeros(EleNo,1);                % Matrix that will hold the mesh elements area
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
        Sta = find(angle(pos_Sta) < ((geo.Qs/3*2-1)*pc*pi/180));
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
        fIron =zeros(IronNo,6);            % Matrix that save flux density of iron element
    end
    
    for kk = 1:1:5
        if jj == 1+(kk-1)*round(nsim/5)               % 5 positions are simulated for losses
            RotPos = exp(j*th_m*pi/180);              % It is a parameter that considers the rotor position
            for i = 1:IronNo
                if (groNoIron(i) == 22)                       % rotor iron elements
                    Pos_Rot = posIron(i)*RotPos;              % get element position
                    %% axis transform
                    fIron(i,kk) = (mo_getb(real(Pos_Rot),imag(Pos_Rot))*[1;j]);
                    Bd = real(fIron(i,kk)) * cosd(th_m) + imag(fIron(i,kk)) * sind(th_m);
                    Bq = imag(fIron(i,kk)) * cosd(th_m) - real(fIron(i,kk)) * sind(th_m);
                    fIron(i,kk) = [Bd,Bq] * [1;j];
                elseif (groNoIron(i) == 12)                   % stator iron elements
                    if ceil(th_m) < 60/p                      % for stator, just 60 ele degree rotation data is needed
                        %%          first one third of stator
                        Pos_Sta = posIron(i);
                        fIron(i,kk) = (mo_getb(real(Pos_Sta),imag(Pos_Sta))*[1;j]);
                        %%          second one third stator iron area
                        Pos_Sta1 = Pos_Sta*exp(j*(geo.Qs/3*2*pc)*pi/180);
                        fIron(i,kk+6) = (mo_getb(real(Pos_Sta1),imag(Pos_Sta1))*[1;j]);
                        %%          last one third stator iron area
                        Pos_Sta2 = Pos_Sta1*exp(j*(geo.Qs/3*2*pc)*pi/180);
                        fIron(i,kk+6*2) = (mo_getb(real(Pos_Sta2),imag(Pos_Sta2))*[1;j]);
                    end
                end
            end
        end
    end
    
    if jj == nsim                                 % Last postion data collection
        RotPos = exp(j*th_m*pi/180);              % It is a parameter that considers the rotor position
        for i = 1:IronNo
            if (groNoIron(i) == 22)                       % rotor iron elements
                Pos_Rot = posIron(i)*RotPos;              % get element position
                %% axis transform
                fIron(i,6) = (mo_getb(real(Pos_Rot),imag(Pos_Rot))*[1;j]);
                Bd = real(fIron(i,6)) * cosd(th_m) + imag(fIron(i,6)) * sind(th_m);
                Bq = imag(fIron(i,6)) * cosd(th_m) - real(fIron(i,6)) * sind(th_m);
                fIron(i,6) = [Bd,Bq] * [1;j];
            elseif (groNoIron(i) == 12)                   % stator iron elements
                if th_m < 60/p                      % for stator, just 60 ele degree rotation data is needed
                    %%          first one third of stator
                    Pos_Sta = posIron(i);
                    fIron(i,6) = (mo_getb(real(Pos_Sta),imag(Pos_Sta))*[1;j]);
                    %%          second one third stator iron area
                    Pos_Sta1 = Pos_Sta*exp(j*(geo.Qs/3*2*pc)*pi/180);
                    fIron(i,6+6) = (mo_getb(real(Pos_Sta1),imag(Pos_Sta1))*[1;j]);
                    %%          last one third stator iron area
                    Pos_Sta2 = Pos_Sta1*exp(j*(geo.Qs/3*2*pc)*pi/180);
                    fIron(i,6+6*2) = (mo_getb(real(Pos_Sta2),imag(Pos_Sta2))*[1;j]);
                end
            end
        end
        fluxdens = [posIron areaIron groNoIron fIron];      % Saving all the info necessary to calculate Iron Losses
    end
    
end