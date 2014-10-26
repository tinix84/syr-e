

function data = buildDefaultRQ(dataSet,bounds)
% s = dataSet.TypeOfRotor;
% data = [0.45 0.22*ones(1,dataSet.NumOfLayers-1) 0.5*ones(1,dataSet.NumOfLayers)];
% if strcmp(s,'Fluid') || strcmp(s,'Seg')
%     data = [data 0.5*ones(1,dataSet.NumOfLayers)];
% end
[bounds, geo, per] = data0(dataSet);

first_index = 1;
if strcmp(geo.RQnames{first_index},'dalpha_1')
    last_index = first_index + geo.nlay-1;
    data = [0.45 0.22*ones(1,dataSet.NumOfLayers-1)];
else
    last_index = 0;
    data = [];
end

first_index = last_index + 1;
if strcmp(geo.RQnames{first_index},'hc')
    last_index = first_index + geo.nlay - 1;
    data = [data 0.5*ones(1,dataSet.NumOfLayers)];
end

first_index = last_index + 1;
if strcmp(geo.RQnames{first_index},'dx')
    % dx per unit
    last_index = first_index + geo.nlay - 1;
    data = [data 0.*ones(1,dataSet.NumOfLayers)];
end

if length(dataSet.RQnames)>(length(data)+1)
    tempRQ = mean(bounds,2);
    for k = length(data)+1:length(dataSet.RQnames)-1
        eval(['data = [data ' num2str(tempRQ(k)) '];'])
    end
end
data = [data 60];

