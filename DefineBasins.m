function Basins = DefineBasins

% Reads IMBIE basin shapes and combines them into larger basins as specified 
% below. The output is a structure with variables Name, ID, X and Y

% read IMBIE data
BasinShape = shaperead([getenv("froot_data"),'/Antarctica_DrainageBasins_IMBIE/Basins_Antarctica_v02_JDR.shp']);

% specify basins to combine
Basins(1).Name = 'AS';
Basins(1).ID = [193 63 145 64 137 136 62 138 192];
Basins(2).Name = 'PIG';
Basins(2).ID = [193 63 145]; % these are from https://nsidc.org/data/NSIDC-0709/versions/2
Basins(3).Name = 'TW';
Basins(3).ID = [64];
Basins(4).Name = 'CR';
Basins(4).ID = [137 136 62 138];
Basins(5).Name = 'DT';
Basins(5).ID = [192];
Basins(6).Name = 'CD';
Basins(6).ID = [137 136 62 138 192];

% combine IMBIE basins as specified by ID vectors
for ii=1:length(Basins)
    for jj=1:length(Basins(ii).ID)
       X = double(BasinShape(Basins(ii).ID(jj)).X);
       Y = double(BasinShape(Basins(ii).ID(jj)).Y);
       C = unique([X(1:end-1)',Y(1:end-1)'],'rows','stable');       
       PolyShapeVectors(jj) = polyshape(C(:,1),C(:,2));
    end
   PolyOut = union(PolyShapeVectors); clear PolyShapeVectors;
   Basins(ii).X = PolyOut.Vertices(:,1);
   Basins(ii).Y = PolyOut.Vertices(:,2);
end