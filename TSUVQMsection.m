function TSUVQMsection(depth,figurefiletag,runID,makemovie)

%% calculates and saves MITgcm diagnostics:
%% 1.  Section diagnostics
%
%       1.1.1 section.xmid
%       1.1.2 section.ymid
%       1.1.4 section.segmentlength
%       1.1.5 section.wct
%       1.1.6 section.wct_belowz
%       1.1.7 section.wct_abovez
%
%       1.2.1 section.monthly.velperp_depthmean
%       1.2.2 section.monthly.velperp_in_depthmean
%       1.2.3 section.monthly.heatflux_depthintegral
%       1.2.4 section.monthly.heatflux_in_depthintegral
%       1.2.5 section.monthly.heatflux_belowz_depthintegral
%       1.2.6 section.monthly.heatflux_belowz_in_depthintegral
%       1.2.7 section.monthly.heatflux_abovez_depthintegral
%       1.2.8 section.monthly.heatflux_abovez_in_depthintegral
%       1.2.9 section.monthly.volumeflux_depthintegral
%       1.2.10 section.monthly.volumeflux_in_depthintegral
%       1.2.11 section.monthly.volumeflux_belowz_depthintegral
%       1.2.12 section.monthly.volumeflux_belowz_in_depthintegral
%       1.2.13 section.monthly.volumeflux_abovez_depthintegral
%       1.2.14 section.monthly.volumeflux_abovez_in_depthintegral
%       1.2.15 section.monthly.TminTf_depthmean
%       1.2.16 section.monthly.TminTf_in_depthmean
%       1.2.17 section.monthly.TminTf_belowz_depthmean
%       1.2.18 section.monthly.TminTf_belowz_in_depthmean
%       1.2.19 section.monthly.TminTf_abovez_depthmean
%       1.2.20 section.monthly.TminTf_abovez_in_depthmean
%       1.2.21 section.monthly.maxoverturning.#basin ( max_z(|\int_B^z \int_xsect u_perp dx dz|) )
%       1.2.22 section.monthly.overturning.#basin ( \int_B^z \int_xsect u_perp dx dz )
%       
%       1.3.1 section.snapshot.xxx
%
%% 2.  Integral diagnostics
%
%% monthly values
%       2.1.1 integral2D.#basin.monthly.melt_integral
%       2.1.2 integral2D.#basin.monthly.ISarea_integral
%       2.1.3 integral2D.#basin.monthly.TfminTf_mean
%       2.1.4 integral2D.#basin.monthly.bgradx_mean
%       2.1.5 integral2D.#basin.monthly.bgrady_mean
% barotropic flow
%       2.1.6 integral2D.#basin.monthly.btUVEL_mean
%       2.1.7 integral2D.#basin.monthly.btVVEL_mean
%       2.1.8 integral2D.#basin.monthly.btVEL_mean
%       2.1.11 integral2D.#basin.monthly.btVEL_max
% ice-ocean boundary flow
%       2.1.12 integral2D.#basin.monthly.blUVEL_mean
%       2.1.13 integral2D.#basin.monthly.blVVEL_mean
%       2.1.14 integral2D.#basin.monthly.blVEL_mean
%       2.1.15 integral2D.#basin.monthly.blVEL_max
%       2.1.16 integral2D.#basin.monthly.blUStar_mean
% 'ambient' flow
%       2.1.17 integral2D.#basin.monthly.aUVEL_mean
%       2.1.18 integral2D.#basin.monthly.aVVEL_mean
%       2.1.19 integral2D.#basin.monthly.aVEL_mean
%       2.1.20 integral2D.#basin.monthly.aVEL_max
% bottom flow
%       2.1.21 integral2D.#basin.monthly.UVELbottom_mean
%       2.1.22 integral2D.#basin.monthly.VVELbottom_mean
% shear in the flow between boundary layer and ambient ocean
%       2.1.23 integral2D.#basin.monthly.UVEL_SHEAR
%       2.1.24 integral2D.#basin.monthly.VVEL_SHEAR
%       2.1.25 integral2D.#basin.monthly.VEL_SHEAR
%
%% snapshots
% barotropic flow
%       2.2.1 integral2D.#basin.snapshot.btUVEL_mean
%       2.2.2 integral2D.#basin.snapshot.btVVEL_mean
%       2.2.3 integral2D.#basin.snapshot.btVEL_mean
%       2.2.4 integral2D.#basin.snapshot.btVEL_max
% ice-ocean boundary flow
%       2.2.5 integral2D.#basin.snapshot.blUVEL_mean
%       2.2.6 integral2D.#basin.snapshot.blVVEL_mean
%       2.2.7 integral2D.#basin.snapshot.blVEL_mean
%       2.2.8 integral2D.#basin.snapshot.blVEL_max
%       2.2.9 integral2D.#basin.snapshot.blUStar_mean
% ambient flow
%       2.2.10 integral2D.#basin.snapshot.aUVEL_mean
%       2.2.11 integral2D.#basin.snapshot.aVVEL_mean
%       2.2.12 integral2D.#basin.snapshot.aVEL_mean
%       2.2.13 integral2D.#basin.snapshot.aVEL_max
% bottom flow
%       2.2.14 integral2D.#basin.snapshot.UVELbottom_mean
%       2.2.15 integral2D.#basin.snapshot.VVELbottom_mean

if nargin<4
    makemovie = 0;
end

if nargin<2
    runID = "PTDC_003";
    figurefiletag = "IceFront_all"; %IceFront_all, moving500mdraft, fixed500mdraft, boundary
    makemovie = 1;
    depth = -400;
end

if exist("HeatVolumeTransport_"+figurefiletag+"_"+runID+".mat")==2
    Data=load("HeatVolumeTransport_"+figurefiletag+"_"+runID+".mat");    
    dn = Data.dn; %dt = Data.dt;
    %start = Data.start;
    MITTime_existing = Data.MITTime;
    MITTime = Data.MITTime;
    %section = Data.section;
    %integral2D = Data.integral2D;
else
    if contains(runID,["805","905","908","815","825"]) % coupling timestep 3 month
        dn = 4; dt = 3; % for all data: dn=1, dt=[1:3]
    elseif contains(runID,["001","002","003","006"]) % coupling timestep 1 month
        dn = 1; dt = 1; % for all data: dn=1, dt=1
    elseif contains(runID,["825_1month"])
        dn = 12; dt = 1; % for all data: dn=1, dt=1
    elseif contains(runID,["000","004"])
        dn = 1; dt = [1:48];
    end
    MITTime_existing = [];
    start = 1;
end

kk=1;%(start-1)*dt(end)+1;
basins = ["AS","PIG","TW","CR","DT"];

frootm = getenv("froot_uamitgcm")+"cases/"+runID;

%interp = dir([frootm,"/Data/GriddedInterpolants_sBh_Bedmap2.mat"]);
%load([interp.folder,"/",interp.name],"FB");
load(getenv("froot_uamitgcm")+"/Ua_InputData/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat","FB");
load(getenv("froot_uamitgcm")+"/UaMITgcm_source/example/PTDC_777/ua_custom/BoundaryCoordinates.mat");

% define some constants
rhoConst = 1024;
si2dbar = 1e-4;
gravity = 9.81;
a0 = -0.0573; a1 = 0; a2 = 0; c0 = 0.0832; b0 = -7.53e-4; % coefficients for local freezing temp of water
%Cp = 3968; % from Whalin et al., 2020, Nature
Cp = 3974; %MITgcm
Lf = 334e3; %J/kg
        
%% read MITgcm output folders
subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = string({subd(isub).name});
nameFolds(ismember(nameFolds,[".",".."])) = [];

%% setup for plotting
folder = "../Figures/"+runID;
if ~exist(folder)
    mkdir(folder);
end

if makemovie
    vidObj = VideoWriter(folder+"/TSUQMsection_bt_"+figurefiletag+"_"+runID);
    %else
    %    vidObj = VideoWriter(folder+"/TSUQMsection_bt_"+figurefiletag+"_"+runID,"MPEG-4");
    %end
    vidObj.FrameRate = 20; vidObj.Quality = 100;
    open(vidObj);

    Hfig=fig("units","inches","width",100*12/72.27,"height",50*12/72.27,"fontsize",12,"font","Helvetica","defaultAxesColorOrder",[0 0 0; 0 0 0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% M A I N    L O O P %%
%%%%%%%%%%%%%%%%%%%%%%%%
startInd = 1;
endInd = min(225*12,numel(nameFolds));
Iinit = 0;

for ii=startInd:dn:endInd

    %% collect some folder and filenames
	MITpath = frootm+"/output/"+nameFolds(ii)+"/MITgcm";
    MITfile = MITpath+"/output.nc";
    Uapath = frootm+"/output/"+nameFolds(ii)+"/Ua";
    if exist(Uapath,'dir')
        Uafiles = dir(Uapath+"/UaDefaultRun*.mat");
        Uafile = Uafiles(1).folder+"/"+Uafiles(1).name;
    else
        Uafile = frootm+"/ua_custom/"+runID+"-RestartFile.mat";
    end
    
    %% read time epoch MITgcm    
    Time = ncread(MITfile,"time");

    %% extract MITgcm data for each timestep and calculate fluxes
    for tt=1:numel(Time)
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %% TIME AND GEOMETRY %%
        %%%%%%%%%%%%%%%%%%%%%%%
        
        attvalue=ncreadatt(MITfile,"time","units");
        if strfind(attvalue,"seconds")
            epoch = erase(attvalue,"seconds since ");
            epochnum = datenum(epoch);
            MITTime(kk) = double(epochnum+Time(tt)/(24*60*60));
        elseif strfind(attvalue,"days")
            epoch = erase(attvalue,"days since ");
            epochnum = datenum(epoch);
            MITTime(kk) = double(epochnum+Time(tt));
        else
            error("I do not recognise the time format in output.nc");
        end       

        Iexisting = find(MITTime(kk)-MITTime_existing==0);

        if isempty(Iexisting)

            %% Define geometric variables
            lon=double(ncread(MITfile,"XC"));
            lat=double(ncread(MITfile,"YC")); 
            z=-double(ncread(MITfile,"Z"));
            if size(lon,2)>1
                lon = lon(:,1);
                lat = lat(1,:);
            end
            nx=numel(lon); ny=numel(lat); nz=numel(z);
            [LAT,LON,Z] = meshgrid(lat,lon,z); dlon=LON(2,1)-LON(1,1); dlat=LAT(1,2)-LAT(1,1);
            [LON2,LAT2] = ndgrid(lon,lat);
            [LAT2_m,LON2_m] = meshgrid(lat,lon);
        
            [~,Idepth] = min(abs(-z-depth));
            
            %% Ice shelf geometry in MITgcm
            [LON_draft,LAT_draft,~,~,~,~,~,~,~,draft] = PlotMeltRates(MITpath,1,[0 0 0 0 0 0 0 1]);
            draft(draft==0)=NaN;
            FbMIT = griddedInterpolant(LON_draft,LAT_draft,draft);
            
            %% define some variables and make plots at the start
            if Iinit==0
                
                % define section of interest      
                if strfind(Uafile,"RestartFile")
                    load(Uafile,"MUA","F","CtrlVarInRestartFile","GF");
                    CtrlVar = CtrlVarInRestartFile;
                    B=F.B;
                else
                    load(Uafile,"MUA","B","CtrlVar","GF");
                end
                CtrlVar.PlotNodes=0; CtrlVar.PlotXYscale=1e3;
                
                G=fig("units","inches","width",80*12/72.27,"fontsize",16,"font","Helvetica","defaultAxesColorOrder",[0 0 0; 0 0 0]);
                subplot("position",[0.1 0.1 0.75 0.85]); hold on;
                %PlotMuaMesh(CtrlVar,MUA);
                PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates/1e3,B);
                caxis([-1500 -500]); colormap("parula");
                PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"m");
                plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");
                plot([lon(1) lon(end) lon(end) lon(1) lon(1)]/1e3,[lat(1) lat(1) lat(end) lat(end) lat(1)]/1e3,"-b");
                xlim([-1710 -1455]); ylim([-710 -225]);
        
                if contains(figurefiletag,"icefront")
                    x_bound=MUA.coordinates(MUA.Boundary.Nodes,1); y_bound=MUA.coordinates(MUA.Boundary.Nodes,2);
                    I_CF = intersect(MUA.Boundary.Nodes,find(GF.node<1));
                    I_CFedges = find(ismember(MUA.Boundary.Edges(:,1),I_CF) | ismember(MUA.Boundary.Edges(:,2),I_CF));
                    xCF = MUA.coordinates(MUA.Boundary.Edges(I_CFedges,:),1); lonsect = reshape(xCF,numel(I_CFedges),2);
                    yCF = MUA.coordinates(MUA.Boundary.Edges(I_CFedges,:),2); latsect = reshape(yCF,numel(I_CFedges),2);
                    lonsect = flip(lonsect,1); lonsect = flip(lonsect,2);
                    latsect = flip(latsect,1); latsect = flip(latsect,2);
                elseif contains(figurefiletag,"draft")
                    draft = FbMIT(LON2,LAT2);
                    draft_tmp = draft; draft_tmp(draft<=depth)=0; draft_tmp(draft>depth)=1;
                    Ccont_s=contour(LON2,LAT2,draft_tmp,[0.99 0.99]);
                    draft_tmp(isnan(draft))=1;
                    Ccont_a=contour(LON2,LAT2,draft_tmp,[0.99 0.99]);
                    [xcont_s,ycont_s] = C2xyz(Ccont_s); lonsect=[]; latsect=[];
                    [xcont_a,ycont_a] = C2xyz(Ccont_a);
                    mm=1; Areas=[];
                    for cc=1:length(xcont_s)
                        xtemp_s = xcont_s{cc}(:); ytemp_s = ycont_s{cc}(:);
                        dends = sqrt((xtemp_s(1)-xtemp_s(end)).^2+(ytemp_s(1)-ytemp_s(end)).^2);
                        if numel(xtemp_s)>20 && dends>5e3
                            lonsect = [lonsect; [xtemp_s(1:end-1) xtemp_s(2:end)]];
                            latsect = [latsect; [ytemp_s(1:end-1) ytemp_s(2:end)]];
                            totdist = [];
                            %% find matching area contour
                            for nn=1:numel(xcont_a)
                                xtemp_a = xcont_a{nn}(:); ytemp_a = ycont_a{nn}(:);
                                [D,~] = pdist2([xtemp_a(:),ytemp_a(:)],[xtemp_s(:) ytemp_s(:)],'euclidean','Smallest',1);
                                totdist(nn) = sum(D,'all');
                            end
                            [~,Iarea] = min(totdist);
                            Areas(mm).lon = xcont_a{Iarea}(:); 
                            Areas(mm).lat = ycont_a{Iarea}(:);
%                             figure(111); hold on;
%                             plot(xtemp_s,ytemp_s,'ok');
%                             plot(DI(mm).lon,DI(mm).lat,'-r'); 
                            mm=mm+1;
                        end  
                    end
                else
                    coordinates_input = ginput(2)*1e3;
                    x1 = coordinates_input(1,1); y1 = coordinates_input(1,2);
                    x2 = coordinates_input(2,1); y2 = coordinates_input(2,2);
                    L_temp = sqrt((x2-x1).^2+(y2-y1).^2); dL_temp = 0.5e3; n = round(L_temp/dL_temp);
                    xnodes = linspace(x1,x2,n+1); ynodes = linspace(y1,y2,n+1); 
                    lonsect = [xnodes(1:end-1)' xnodes(2:end)']; latsect = [ynodes(1:end-1)' ynodes(2:end)']; 
                end
                
                dL = sqrt((lonsect(:,2)-lonsect(:,1)).^2+(latsect(:,2)-latsect(:,1)).^2);
                d = cumsum(dL);
                epar = [(lonsect(:,2)-lonsect(:,1))./dL,(latsect(:,2)-latsect(:,1))./dL];
                eperp = [-epar(:,2),epar(:,1)];      
                lonmid = 0.5*(lonsect(:,1)+lonsect(:,2));
                latmid = 0.5*(latsect(:,1)+latsect(:,2));

                %figure(111); hold on;
                %scatter(lonmid,latmid,[],d,"filled");
        
                LONSECT = repmat(lonmid,1,numel(z)); 
                LATSECT = repmat(latmid,1,numel(z)); 
                ZSECT = repmat(z(:)',size(lonsect,1),1);
                [Dgrid,Zgrid] = ndgrid(d/1e3,-z);
                
                [inUa,onUa] = inpoly([lonsect(:,1) latsect(:,1)],[MUA.Boundary.x MUA.Boundary.y]);
        
                plot(lonsect'/1e3,latsect'/1e3,"-k");
                %quiver(lonmid/1e3,latmid/1e3,eperp(:,1)*5,eperp(:,2)*5,"linewidth",1,"color","k","autoscale","off");
                %quiver([x1+x2]/2e3,[y1+y2]/2e3,epar(1)*10,epar(2)*10,"linewidth",2,"color","b");
                pos = get(G,"Position");
                set(G,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
                fname = folder+"/Mesh_"+figurefiletag+"_"+runID;
                %print(G,fname,"-dpng","-r400");

                Iinit = 1;
                
            end
            
            % if draft changes: need to recompute gates every step
            if contains(figurefiletag,"moving")
                draft = FbMIT(LON2,LAT2);
                draft_tmp = draft; draft_tmp(draft<=depth)=0; draft_tmp(draft>depth)=1;
                Ccont_s=contour(LON2,LAT2,draft_tmp,[0.99 0.99]);
                draft_tmp(isnan(draft))=1;
                Ccont_a=contour(LON2,LAT2,draft_tmp,[0.99 0.99]);
                [xcont_s,ycont_s] = C2xyz(Ccont_s); lonsect=[]; latsect=[];
                [xcont_a,ycont_a] = C2xyz(Ccont_a);
                mm=1; Areas=[];
                for cc=1:length(xcont_s)
                    xtemp_s = xcont_s{cc}(:); ytemp_s = ycont_s{cc}(:);
                    dends = sqrt((xtemp_s(1)-xtemp_s(end)).^2+(ytemp_s(1)-ytemp_s(end)).^2);
                    if numel(xtemp_s)>20 && dends>5e3
                        lonsect = [lonsect; [xtemp_s(1:end-1) xtemp_s(2:end)]];
                        latsect = [latsect; [ytemp_s(1:end-1) ytemp_s(2:end)]];
                        totdist = [];
                        %% find matching area contour
                        for nn=1:numel(xcont_a)
                            xtemp_a = xcont_a{nn}(:); ytemp_a = ycont_a{nn}(:);
                            [D,~] = pdist2([xtemp_a(:),ytemp_a(:)],[xtemp_s(:) ytemp_s(:)],'euclidean','Smallest',1);
                            totdist(nn) = sum(D,'all');
                        end
                        [~,Iarea] = min(totdist);
                        Areas(mm).lon = xcont_a{Iarea}(:); 
                        Areas(mm).lat = ycont_a{Iarea}(:);
                        %figure(111); hold on;
                        %plot(xtemp_s,ytemp_s,'ok');
                        %plot(DI(mm).lon,DI(mm).lat,'-r'); 
                        mm=mm+1;
                    end  
                end
                dL = sqrt((lonsect(:,2)-lonsect(:,1)).^2+(latsect(:,2)-latsect(:,1)).^2);
                d = cumsum(dL);
                epar = [(lonsect(:,2)-lonsect(:,1))./dL,(latsect(:,2)-latsect(:,1))./dL];
                eperp = [-epar(:,2),epar(:,1)];      
                lonmid = 0.5*(lonsect(:,1)+lonsect(:,2));
                latmid = 0.5*(latsect(:,1)+latsect(:,2));
                %figure; scatter(lonmid,latmid,[],d,"filled");
                
                LONSECT = repmat(lonmid,1,numel(z)); 
                LATSECT = repmat(latmid,1,numel(z)); 
                ZSECT = repmat(z(:)',size(lonsect,1),1);
                [Dgrid,Zgrid] = ndgrid(d/1e3,-z);
                
                [inUa,onUa] = inpoly([lonsect(:,1) latsect(:,1)],[MUA.Boundary.x MUA.Boundary.y]);
                
        %         if mod(ii,dn)==0
        %             figure(G); plot(lonsect"/1e3,latsect"/1e3,"-k");
        %         end
            end
            
            %% extract geometry and melt data along section for plotting
            if makemovie
                if strfind(Uafile,"RestartFile")
                    load(Uafile,"MUA","F","CtrlVarInRestartFile","GF");
                    CtrlVar = CtrlVarInRestartFile;
                    b = F.b;
                    s = F.s;
                    B = F.B;
                    ab = F.ab;
                    ub = F.ub;
                    vb = F.vb;
                else
                    load(Uafile,"MUA","b","s","B","CtrlVar","GF","ab","ub","vb");  
                end
                
                x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);
                Fb = scatteredInterpolant(x,y,b);
                b(GF.node==1)=NaN;
                %FbIS = scatteredInterpolant(x,y,b); 
                Fs = scatteredInterpolant(x,y,s);
                Fab =  scatteredInterpolant(x,y,ab);
                Fspeed = scatteredInterpolant(x,y,sqrt(ub.^2+vb.^2));
                bsect = Fb(lonmid,latmid); bsect(inUa==0)=NaN;
                ssect = Fs(lonmid,latmid); ssect(inUa==0)=NaN;
                usect = Fspeed(lonmid,latmid); usect(inUa==0)=NaN;
                Bsect = FB(lonmid,latmid);
                absect = Fab(lonmid,latmid); absect(inUa==0)=NaN;
                if ii==start
                   bsect_init = bsect;
                   ssect_init = ssect;
                   usect_init = usect;
                   d_init = d;
                end
            end

            drF = ncread(MITfile,"drF");
            drF_full = permute(repmat(drF,1,nx,ny),[2 3 1]); 
            
            hFacS_full = ncread(MITfile,"hFacS"); 
            dyG = ncread(MITfile,"dyG");
            dyG_full = repmat(dyG,1,1,nz);
            
            hFacW_full = ncread(MITfile,"hFacW");
            dxG = ncread(MITfile,"dxG"); 
            dxG_full = repmat(dxG,1,1,nz);
            
            hFacC_full = ncread(MITfile,"hFacC");
            
            depthC = cumsum(ncread(MITfile,"drC"));
            p_depth = depthC(1:end-1)*gravity*rhoConst*si2dbar;
            p_full= permute(repmat(p_depth,1,nx,ny),[2 3 1]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            %% SECTION DIAGNOSTICS %%
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            % load 3D fields for U, V, T and S (monthly averages)
            UVEL=ncread(MITfile,"UVEL",[1,1,1,tt],[Inf,Inf,Inf,1]);
            VVEL=ncread(MITfile,"VVEL",[1,1,1,tt],[Inf,Inf,Inf,1]);
            THETA=ncread(MITfile,"THETA",[1,1,1,tt],[Inf,Inf,Inf,1]);        
            SALT=ncread(MITfile,"SALT",[1,1,1,tt],[Inf,Inf,Inf,1]);
            THETA = gsw_pt_from_t(SALT,THETA,0*THETA,p_full); %turn potential temperature into in-situ temperature relative to the surface
            
            Rho =  densmdjwf(SALT,THETA,p_full);
            %Rho_temp = Rho; Rho_temp(Rho_temp==0)=NaN;
            %Rho_depthav = mean(Rho_temp,3,"omitnan");
            Tfreeze = a0*SALT+c0+b0*p_full;
            DelTHETA = THETA-Tfreeze;
            
            H = CalcWCTMIT(MITfile,tt);
            [~,~,zbottom,~,~,~,~] = BottomProperties(MITfile,tt);
            
            % calculate velocity components perpendicular to section
            UVEL_sect = interp3(LAT,LON,Z,UVEL,LATSECT,LONSECT,ZSECT);
            VVEL_sect = interp3(LAT,LON,Z,VVEL,LATSECT,LONSECT,ZSECT);
            VEL_perpsect = UVEL_sect.*repmat(eperp(:,1),1,length(z)) + VVEL_sect.*repmat(eperp(:,2),1,length(z));
            Iout = find(VEL_perpsect<0); % find outflow velocities
            VEL_perpsect(VEL_perpsect==0)=NaN; % set zero velocities to NaN for averaging
            VEL_perpsect_in = VEL_perpsect; VEL_perpsect_in(Iout)=NaN;
            
            THETA_sect = interp3(LAT,LON,Z,THETA,LATSECT,LONSECT,ZSECT);   
            THETA_sect(isnan(VEL_perpsect))=NaN;
            THETA_sect_in = THETA_sect; THETA_sect_in(Iout)=NaN;
            SALT_sect = interp3(LAT,LON,Z,SALT,LATSECT,LONSECT,ZSECT);
            SALT_sect(isnan(VEL_perpsect))=NaN;
            SALT_sect_in = SALT_sect; SALT_sect_in(Iout)=NaN;
                
            UVEL_Wfaces = UVEL.*drF_full.*hFacW_full; %in m2/s
            VVEL_Sfaces = VVEL.*drF_full.*hFacS_full; %in m2/s
            
            HeatFlux_Wfaces = Cp*Rho.*UVEL_Wfaces.*DelTHETA/1e12; % in J/(m*s)
            HeatFlux_Sfaces = Cp*Rho.*VVEL_Sfaces.*DelTHETA/1e12;   
            HeatFlux_Wfaces_sect = interp3(LAT,LON,Z,HeatFlux_Wfaces,LATSECT,LONSECT,ZSECT);
            HeatFlux_Sfaces_sect = interp3(LAT,LON,Z,HeatFlux_Sfaces,LATSECT,LONSECT,ZSECT);
            HeatFlux_sect = (HeatFlux_Wfaces_sect.*repmat(eperp(:,1),1,length(z)) + ...
                HeatFlux_Sfaces_sect.*repmat(eperp(:,2),1,length(z))).*repmat(dL,1,length(z)); % in J/s     
            HeatFlux_sect_in = HeatFlux_sect; HeatFlux_sect_in(Iout)=0;
            
            VolumeFlux_Wfaces_sect = interp3(LAT,LON,Z,UVEL_Wfaces,LATSECT,LONSECT,ZSECT);
            VolumeFlux_Sfaces_sect = interp3(LAT,LON,Z,VVEL_Sfaces,LATSECT,LONSECT,ZSECT);
            VolumeFlux_sect = (VolumeFlux_Wfaces_sect.*repmat(eperp(:,1),1,length(z)) + ...
                VolumeFlux_Sfaces_sect.*repmat(eperp(:,2),1,length(z))).*repmat(dL,1,length(z)); % in m3/s
            VolumeFlux_sect_in = VolumeFlux_sect; VolumeFlux_sect_in(Iout)=0;
            
            DelTheta_sect = interp3(LAT,LON,Z,DelTHETA,LATSECT,LONSECT,ZSECT);
            DelTheta_sect(isnan(VEL_perpsect))=NaN;
            DelTheta_sect_in = DelTheta_sect; DelTheta_sect_in(Iout)=NaN;
    
            % assign variables to structure
            % ETAN corrections are ignored for now as they are small
            section(kk).xmid = lonmid;
            section(kk).ymid = latmid;
            section(kk).segmentlength = dL;
            section(kk).wct = interp2(LAT2_m,LON2_m,H,latmid,lonmid);
            section(kk).wct_belowz = max([0*lonmid(:),depth(:)-interp2(LAT2_m,LON2_m,zbottom,latmid(:),lonmid(:))],[],2);
            section(kk).wct_abovez = max([0*lonmid(:),(interp2(LAT2,LON2,zbottom,latmid(:),lonmid(:))+section(kk).wct(:))-depth(:)],[],2);
            
            section(kk).monthly.velperp_depthmean = mean(VEL_perpsect,2,"omitnan");
            section(kk).monthly.velperp_in_depthmean = mean(VEL_perpsect_in,2,"omitnan");
            section(kk).monthly.heatflux_depthintegral = sum(HeatFlux_sect,2,"omitnan");
            section(kk).monthly.heatflux_in_depthintegral = sum(HeatFlux_sect_in,2,"omitnan");
            section(kk).monthly.heatflux_belowz_depthintegral = sum(HeatFlux_sect(:,Idepth:end),2,"omitnan");
            section(kk).monthly.heatflux_belowz_in_depthintegral = sum(HeatFlux_sect_in(:,Idepth:end),2,"omitnan");
            section(kk).monthly.heatflux_abovez_depthintegral = sum(HeatFlux_sect(:,1:Idepth-1),2,"omitnan");
            section(kk).monthly.heatflux_abovez_in_depthintegral = sum(HeatFlux_sect_in(:,1:Idepth-1),2,"omitnan");
            
            section(kk).monthly.volumeflux_depthintegral = sum(VolumeFlux_sect,2,"omitnan");
            section(kk).monthly.volumeflux_in_depthintegral = sum(VolumeFlux_sect_in,2,"omitnan");
            section(kk).monthly.volumeflux_belowz_depthintegral = sum(VolumeFlux_sect(:,Idepth:end),2,"omitnan");
            section(kk).monthly.volumeflux_belowz_in_depthintegral = sum(VolumeFlux_sect_in(:,Idepth:end),2,"omitnan");
            section(kk).monthly.volumeflux_abovez_depthintegral = sum(VolumeFlux_sect(:,1:Idepth-1),2,"omitnan");
            section(kk).monthly.volumeflux_abovez_in_depthintegral = sum(VolumeFlux_sect_in(:,1:Idepth-1),2,"omitnan");
            
            section(kk).monthly.TminTf_depthmean = mean(DelTheta_sect,2,"omitnan");
            section(kk).monthly.TminTf_in_depthmean = mean(DelTheta_sect_in,2,"omitnan");
            section(kk).monthly.TminTf_belowz_depthmean = mean(DelTheta_sect(:,Idepth:end),2,"omitnan");
            section(kk).monthly.TminTf_belowz_in_depthmean = mean(DelTheta_sect_in(:,Idepth:end),2,"omitnan");
            section(kk).monthly.TminTf_abovez_depthmean = mean(DelTheta_sect(:,1:Idepth-1),2,"omitnan");
            section(kk).monthly.TminTf_abovez_in_depthmean = mean(DelTheta_sect_in(:,1:Idepth-1),2,"omitnan");
    
            section(kk).monthly.T_depthmean = mean(THETA_sect,2,"omitnan");
            section(kk).monthly.T_in_depthmean = mean(THETA_sect_in,2,"omitnan");
            section(kk).monthly.T_belowz_depthmean = mean(THETA_sect(:,Idepth:end),2,"omitnan");
            section(kk).monthly.T_belowz_in_depthmean = mean(THETA_sect_in(:,Idepth:end),2,"omitnan");
            section(kk).monthly.T_abovez_depthmean = mean(THETA_sect(:,1:Idepth-1),2,"omitnan");
            section(kk).monthly.T_abovez_in_depthmean = mean(THETA_sect_in(:,1:Idepth-1),2,"omitnan");
    
            section(kk).monthly.S_depthmean = mean(SALT_sect,2,"omitnan");
            section(kk).monthly.S_in_depthmean = mean(SALT_sect_in,2,"omitnan");
            section(kk).monthly.S_belowz_depthmean = mean(SALT_sect(:,Idepth:end),2,"omitnan");
            section(kk).monthly.S_belowz_in_depthmean = mean(SALT_sect_in(:,Idepth:end),2,"omitnan");
            section(kk).monthly.S_abovez_depthmean = mean(SALT_sect(:,1:Idepth-1),2,"omitnan");
            section(kk).monthly.S_abovez_in_depthmean = mean(SALT_sect_in(:,1:Idepth-1),2,"omitnan");
    
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %% INTEGRAL DIAGNOSTICS %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %% load monthly data
            [~,~,MeltMIT,UVEL_bl,VVEL_bl,UStar,~,~,TfminT,Draft] = PlotMeltRates(MITpath,tt,[1 1 1 1 0 0 1 1]);
            VEL_bl = sqrt(UVEL_bl.^2+VVEL_bl.^2);
            [~,~,~,~,~,~,~,UVEL_bt,VVEL_bt,~,~,~,~,~] = CalcVorticity(MITfile,tt);
            VEL_bt = sqrt(UVEL_bt.^2+VVEL_bt.^2);
            UVEL_tmp = UVEL; VVEL_tmp = VVEL; UVEL_tmp(UVEL_tmp==0)=NaN; VVEL_tmp(VVEL_tmp==0)=NaN;
            dUVEL = abs(UVEL_bl-UVEL_bt);
            dVVEL = abs(VVEL_bl-VVEL_bt);
            
            %dVEL = abs(VEL_bl-VEL_bt);
            %UVELratio = dUVEL./abs(UVEL_bl);
            %VVELratio = dVVEL./abs(VVEL_bl);
            %VELratio = dVEL./abs(VEL_bl);

            %VELratio = dUVEL./UVEL_bl;

    %         UVEL_bc = UVEL - UVEL_bt; VVEL_bc = VVEL - VVEL_bt;
    %         VEL_bc = sqrt(UVEL_bc.^2+VVEL_bc.^2);
            H_bl = min(0*H+z(2)-z(1),H); 
            UVEL_a = (H.*UVEL_bt - H_bl.*UVEL_bl)./(H-H_bl); UVEL_a(abs(UVEL_a)==Inf)=NaN;
            VVEL_a = (H.*VVEL_bt - H_bl.*VVEL_bl)./(H-H_bl); VVEL_a(abs(VVEL_a)==Inf)=NaN;
            VEL_a = sqrt(UVEL_a.^2+VVEL_a.^2);

            Draft(Draft>=0) = NaN;
            [gradby,gradbx] = gradient(Draft,dlat,dlon);

            VELratio_gradb = (dUVEL.*abs(gradbx) + dVVEL.*abs(gradby))./VEL_bl;
            UVELratio = dUVEL./VEL_bl;
            VVELratio = dVVEL./VEL_bl;
            gradbx_abs = abs(gradbx);
            gradby_abs = abs(gradby);

            E0 = 3.6e-2;
            sqCd = sqrt(2.5e-3); %non-dim
            GammT = 0.02;

            epsT = VELratio_gradb*E0./(sqCd*GammT+VELratio_gradb*E0);
            
            %% load snapshot data
            nstr = strlength(MITpath);
            moy = str2double(extractBetween(MITpath,nstr-8,nstr-7));
            if tt==numel(Time) && moy==1
                [UVEL_bt_snapshot,VVEL_bt_snapshot,~,~,UVEL_bottom_snapshot,VVEL_bottom_snapshot] = ...
                        CalcSnapshotVelocityComponents(MITpath);
                UVEL_bt_snapshot = reshape(UVEL_bt_snapshot,nx,ny);
                VVEL_bt_snapshot = reshape(VVEL_bt_snapshot,nx,ny);
                VEL_bt_snapshot = sqrt(UVEL_bt_snapshot.^2+VVEL_bt_snapshot.^2);
                [~,~,~,~,~,UVEL_bl_snapshot,drTop_u,VVEL_bl_snapshot,drTop_v,UStar_snapshot,...
                    UVEL_a_snapshot,drRem_u,VVEL_a_snapshot,drRem_v] = CalcSnapshotMITBoundaryLayerProperties(MITpath); 
                VEL_bl_snapshot = sqrt(UVEL_bl_snapshot.^2+VVEL_bl_snapshot.^2);
                %UVEL_a_snapshot = (H.*UVEL_bt_snapshot - H_bl.*UVEL_bl_snapshot)./(H-H_bl); UVEL_a_snapshot(abs(UVEL_a_snapshot)==Inf)=NaN;
                %VVEL_a_snapshot = (H.*VVEL_bt_snapshot - H_bl.*VVEL_bl_snapshot)./(H-H_bl); VVEL_a_snapshot(abs(VVEL_a_snapshot)==Inf)=NaN;
                VEL_a_snapshot = sqrt(UVEL_a_snapshot.^2+VVEL_a_snapshot.^2);
            else
                VEL_bt_snapshot = NaN*MeltMIT;
                UVEL_bt_snapshot = NaN*MeltMIT;
                VVEL_bt_snapshot = NaN*MeltMIT;
                VEL_bl_snapshot = NaN*MeltMIT;
                UVEL_bl_snapshot = NaN*MeltMIT;
                VVEL_bl_snapshot = NaN*MeltMIT;
                UStar_snapshot = NaN*MeltMIT;
                UVEL_a_snapshot = NaN*MeltMIT;
                VVEL_a_snapshot = NaN*MeltMIT;
                VEL_a_snapshot = NaN*MeltMIT;
    %            UVEL_bottom_snapshot = NaN*MeltMIT;
    %            VVEL_bottom_snapshot = NaN*MeltMIT;
            end
            
            if contains(figurefiletag,"IceFront")
                %MeltMIT [kg/s/m2] * Lf [J/kg] * dx [m] * dy [m] / 1e12
                Imelt = find(MeltMIT~=0);
            elseif contains(figurefiletag,"draft")
                Imelt = [];
                for nn=1:numel(Areas)
                    I_tmp = find(inpoly([LON2(:),LAT2(:)],[Areas(nn).lon,Areas(nn).lat]));
                    Imelt = [Imelt; I_tmp];
                end
                Imelt = unique(Imelt);
            else
                Imelt = [];
            end
            
            if ~isempty(Imelt)
                %% define basins
                % AS
                basin(1).name = "AS";
                basin(1).Imelt = Imelt;
                basin(1).Isect = [1:numel(section(kk).xmid)];

                % PIG
                basin(2).name = "PIG";
                %xmin = -1699e3; xmax = -1500e3; ymin = -375e3; ymax = -220e3;
                %I = find(LON2(Imelt)>xmin & LON2(Imelt)<xmax & LAT2(Imelt)>ymin & LAT2(Imelt)<ymax); 
                %basin(2).Imelt = Imelt(I);
                %basin(2).Isect = find(section(kk).xmid(:)>xmin & section(kk).xmid(:)<xmax & section(kk).ymid(:)>ymin & section(kk).ymid(:)<ymax);
                xpoly = [-1699 -1699 -1500 -1500 -1550]*1e3;
                ypoly = [-375 -220 -220 -300 -375]*1e3;
                I = find(inpoly([LON2(Imelt(:)),LAT2(Imelt(:))],[xpoly(:) ypoly(:)]));
                basin(2).Imelt = Imelt(I); 
                basin(2).Isect = find(inpoly([section(kk).xmid(:) section(kk).ymid(:)],[xpoly(:) ypoly(:)]));
                
                % Thwaites
                basin(3).name = "TW";
                %xmin = -1620e3; xmax = -1500e3; ymin = -520e3; ymax = -375e3;
                %I = find(LON2(Imelt)>xmin & LON2(Imelt)<xmax & LAT2(Imelt)>ymin & LAT2(Imelt)<ymax);
                %basin(3).Imelt = Imelt(I);
                %basin(3).Isect = find(section(kk).xmid(:)>xmin & section(kk).xmid(:)<xmax & section(kk).ymid(:)>ymin & section(kk).ymid(:)<ymax);
                xpoly = [-1620 -1550 -1500 -1450 -1450 -1620]*1e3;
                ypoly = [-375 -375 -300 -300 -520 -520]*1e3;
                I = find(inpoly([LON2(Imelt(:)),LAT2(Imelt(:))],[xpoly(:) ypoly(:)]));
                basin(3).Imelt = Imelt(I); 
                basin(3).Isect = find(inpoly([section(kk).xmid(:) section(kk).ymid(:)],[xpoly(:) ypoly(:)]));
                
                % Crosson
                basin(4).name = "CR";
                xpoly = [-1610 -1485 -1450 -1450 -1610]*1e3;
                ypoly = [-580 -657 -657 -520 -520]*1e3; 
                I = find(inpoly([LON2(Imelt(:)),LAT2(Imelt(:))],[xpoly(:) ypoly(:)]));
                basin(4).Imelt = Imelt(I);
                basin(4).Isect = find(inpoly([section(kk).xmid(:) section(kk).ymid(:)],[xpoly(:) ypoly(:)]));

                % Dotson
                basin(5).name = "DT";
                xpoly = [-1610 -1485 -1450 -1450 -1625]*1e3;
                ypoly = [-580 -657 -657 -700 -700]*1e3; 
                I = find(inpoly([LON2(Imelt(:)),LAT2(Imelt(:))],[xpoly(:) ypoly(:)]));
                basin(5).Imelt = Imelt(I); 
                basin(5).Isect = find(inpoly([section(kk).xmid(:) section(kk).ymid(:)],[xpoly(:) ypoly(:)]));
                
                CM=[1 0 0;0 1 0;0 0 1;1 0 1;0 1 1];
                
                for bb=1:length(basins)

                    BasinVarname=basin(bb).name;
                    Imelt = basin(bb).Imelt;
                    Isect = basin(bb).Isect;
                    %CM=[1 0 0;0 1 0;0 0 1;1 0 1;0 1 1];
                    %figure(G); plot(LON2(Imelt)/1e3,LAT2(Imelt)/1e3,".","color",CM(bb,:));
    
                    %%% MONTHLY AVERAGES %%%
                    integral2D.(BasinVarname).monthly.draft_mean(kk) = mean(Draft(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.melt_integral(kk) =  sum(-MeltMIT(Imelt),"all","omitnan")*(365.25*24*60*60)*dlon*dlat/1e12; % in Gt/yr 
                    integral2D.(BasinVarname).monthly.ISarea_integral(kk) = numel(Imelt)*dlon*dlat;
                    integral2D.(BasinVarname).monthly.bgradx_mean(kk) = mean(gradbx(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.bgrady_mean(kk) = mean(gradby(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.TminTf_mean(kk) = mean(TfminT(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.btVEL_mean(kk) = mean(VEL_bt(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.btUVEL_mean(kk) = mean(UVEL_bt(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.btVVEL_mean(kk) = mean(VVEL_bt(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.blVEL_mean(kk) = mean(VEL_bl(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.blUVEL_mean(kk) = mean(UVEL_bl(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.blVVEL_mean(kk) = mean(VVEL_bl(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.blUStar_mean(kk) = mean(UStar(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.aVEL_mean(kk) = mean(VEL_a(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.aUVEL_mean(kk) = mean(UVEL_a(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.aVVEL_mean(kk) = mean(VVEL_a(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).monthly.VELratio_gradb(kk) = mean(VELratio_gradb(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.UVELratio(kk) = mean(UVELratio(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.VVELratio(kk) = mean(VVELratio(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.gradbx_abs(kk) = mean(gradbx_abs(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.gradby_abs(kk) = mean(gradby_abs(Imelt),"all","omitnan"); 
                    integral2D.(BasinVarname).monthly.epsT(kk) = mean(epsT(Imelt),"all","omitnan"); 

                    if ~isempty(Imelt)
                        integral2D.(BasinVarname).monthly.btVEL_max(kk) = max(VEL_bt(Imelt),[],"all","omitnan");
                        integral2D.(BasinVarname).monthly.blVEL_max(kk) = max(VEL_bl(Imelt),[],"all","omitnan");
                        integral2D.(BasinVarname).monthly.aVEL_max(kk) = max(VEL_a(Imelt),[],"all","omitnan"); 
                    else
                        integral2D.(BasinVarname).monthly.btVEL_max(kk) = NaN;
                        integral2D.(BasinVarname).monthly.blVEL_max(kk) = NaN;
                        integral2D.(BasinVarname).monthly.aVEL_max(kk) = NaN;
                    end

                    
        %            integral2D(kk).(BasinVarname).monthly.UVELbottom_mean = mean(UVEL_bottom(Imelt),"all","omitnan");
        %            integral2D(kk).(BasinVarname).monthly.VVELbottom_mean = mean(VVEL_bottom(Imelt),"all","omitnan");
    
                    %%% SNAPSHOT %%%
                    integral2D.(BasinVarname).snapshot.btVEL_mean(kk) = mean(VEL_bt_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.btUVEL_mean(kk) = mean(UVEL_bt_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.btVVEL_mean(kk) = mean(VVEL_bt_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.blVEL_mean(kk) = mean(VEL_bl_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.blUVEL_mean(kk) = mean(UVEL_bl_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.blVVEL_mean(kk) = mean(VVEL_bl_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.blUStar_mean(kk) = mean(UStar_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.aVEL_mean(kk) = mean(VEL_a_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.aUVEL_mean(kk) = mean(UVEL_a_snapshot(Imelt),"all","omitnan");
                    integral2D.(BasinVarname).snapshot.aVVEL_mean(kk) = mean(VVEL_a_snapshot(Imelt),"all","omitnan");

                    if ~isempty(Imelt)
                        integral2D.(BasinVarname).snapshot.btVEL_max(kk) = max(VEL_bt_snapshot(Imelt),[],"all","omitnan"); 
                        integral2D.(BasinVarname).snapshot.blVEL_max(kk) = max(VEL_bl_snapshot(Imelt),[],"all","omitnan");
                        integral2D.(BasinVarname).snapshot.aVEL_max(kk) = max(VEL_a_snapshot(Imelt),[],"all","omitnan");
                    else
                        integral2D.(BasinVarname).snapshot.btVEL_max(kk) = NaN; 
                        integral2D.(BasinVarname).snapshot.blVEL_max(kk) = NaN;
                        integral2D.(BasinVarname).snapshot.aVEL_max(kk) = NaN;
                    end
                    
        %            integral2D(kk).(BasinVarname).snapshot.UVELbottom_mean = mean(UVEL_bottom_snapshot(Imelt),"all","omitnan");
        %            integral2D(kk).(BasinVarname).snapshot.VVELbottom_mean = mean(VVEL_bottom_snapshot(Imelt),"all","omitnan");
    
                    %%% SECTION OVERTURNING %%%
                    section(kk).monthly.overturning.(BasinVarname)=sum(VolumeFlux_sect(Isect,:),1);
                    section(kk).monthly.maxoverturning.(BasinVarname)=max(abs(cumsum(sum(VolumeFlux_sect(Isect,:),1))));
                    
                end


            end
            
            if contains(figurefiletag,"draft")
                for nn=1:numel(Areas)
                    DI(kk).Areas(nn).lon=Areas(nn).lon;
                    DI(kk).Areas(nn).lat=Areas(nn).lat;
                end
            else
                DI(kk)=NaN;
            end

            titletime{kk} =  datestr(MITTime(kk),"dd-mm-yyyy");
            fprintf("New: %s\n",titletime{kk});

        else

            Iexisting = Iexisting(1);
            section(kk) = Data.section(Iexisting);
            DI(kk) = Data.DI(Iexisting);

            for bb=1:length(basins)
                    
                BasinVarname=basins{bb};

                %%% MONTHLY AVERAGES %%%
                integral2D.(BasinVarname).monthly.draft_mean(kk) = Data.integral2D.(BasinVarname).monthly.draft_mean(Iexisting);
                integral2D.(BasinVarname).monthly.melt_integral(kk) = Data.integral2D.(BasinVarname).monthly.melt_integral(Iexisting);
                integral2D.(BasinVarname).monthly.ISarea_integral(kk) = Data.integral2D.(BasinVarname).monthly.ISarea_integral(Iexisting);
                integral2D.(BasinVarname).monthly.bgradx_mean(kk) =  Data.integral2D.(BasinVarname).monthly.bgradx_mean(Iexisting);
                integral2D.(BasinVarname).monthly.bgrady_mean(kk) = Data.integral2D.(BasinVarname).monthly.bgrady_mean(Iexisting);
                integral2D.(BasinVarname).monthly.TminTf_mean(kk) = Data.integral2D.(BasinVarname).monthly.TminTf_mean(Iexisting);
                integral2D.(BasinVarname).monthly.btVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.btVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.btVEL_max(kk) = Data.integral2D.(BasinVarname).monthly.btVEL_max(Iexisting);
                integral2D.(BasinVarname).monthly.btUVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.btUVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.btVVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.btVVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.blVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.blVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.blVEL_max(kk) = Data.integral2D.(BasinVarname).monthly.blVEL_max(Iexisting);
                integral2D.(BasinVarname).monthly.blUVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.blUVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.blVVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.blVVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.blUStar_mean(kk) = Data.integral2D.(BasinVarname).monthly.blUStar_mean(Iexisting);
                integral2D.(BasinVarname).monthly.aVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.aVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.aVEL_max(kk) = Data.integral2D.(BasinVarname).monthly.aVEL_max(Iexisting);
                integral2D.(BasinVarname).monthly.aUVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.aUVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.aVVEL_mean(kk) = Data.integral2D.(BasinVarname).monthly.aVVEL_mean(Iexisting);
                integral2D.(BasinVarname).monthly.VELratio_gradb(kk) = Data.integral2D.(BasinVarname).monthly.VELratio_gradb(Iexisting);
                integral2D.(BasinVarname).monthly.UVELratio(kk) = Data.integral2D.(BasinVarname).monthly.UVELratio(Iexisting);
                integral2D.(BasinVarname).monthly.VVELratio(kk) = Data.integral2D.(BasinVarname).monthly.VVELratio(Iexisting);
                integral2D.(BasinVarname).monthly.gradbx_abs(kk) = Data.integral2D.(BasinVarname).monthly.gradbx_abs(Iexisting);
                integral2D.(BasinVarname).monthly.gradby_abs(kk) = Data.integral2D.(BasinVarname).monthly.gradby_abs(Iexisting);
                integral2D.(BasinVarname).monthly.epsT(kk) = Data.integral2D.(BasinVarname).monthly.epsT(Iexisting);
    %      
                %%% SNAPSHOT %%%
                integral2D.(BasinVarname).snapshot.btVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.btVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.btVEL_max(kk) = Data.integral2D.(BasinVarname).snapshot.btVEL_max(Iexisting);
                integral2D.(BasinVarname).snapshot.btUVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.btUVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.btVVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.btVVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.blVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.blVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.blVEL_max(kk) = Data.integral2D.(BasinVarname).snapshot.blVEL_max(Iexisting);
                integral2D.(BasinVarname).snapshot.blUVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.blUVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.blVVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.blVVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.blUStar_mean(kk) = Data.integral2D.(BasinVarname).snapshot.blUStar_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.aVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.aVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.aVEL_max(kk) = Data.integral2D.(BasinVarname).snapshot.aVEL_max(Iexisting);
                integral2D.(BasinVarname).snapshot.aUVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.aUVEL_mean(Iexisting);
                integral2D.(BasinVarname).snapshot.aVVEL_mean(kk) = Data.integral2D.(BasinVarname).snapshot.aVVEL_mean(Iexisting);

            end

            titletime{kk} = datestr(MITTime_existing(Iexisting),"dd-mm-yyyy");
            fprintf("Existing: %s\n",titletime{kk});

        end

        if makemovie
                %% add new frame to video
                figure(Hfig);
    
                subplot("position",[0.06 0.75 0.3 0.2]); hold on;
    
                abmelt = absect; abmelt(abmelt>0)=NaN;
                abfreeze = absect; abfreeze(abfreeze<0)=NaN;
                plot(d/1e3,abmelt,"-r");
                plot(d/1e3,abfreeze,"-b");
                xlim([0 max(d)/1e3-1]); xticklabels({});
                ylim([-150 0]); ylabel("Melt [m/yr]");
                grid on; box on;
    
                subplot("position",[0.06 0.1 0.3 0.595]); hold on;
    
                pcolor(Dgrid,Zgrid,THETA_sect); shading flat;
                SALT_sect(SALT_sect>34.65 | SALT_sect<33.95)=NaN;
                [C,h]=contour(Dgrid,Zgrid,SALT_sect,[34 34.3 34.6],"color","k","linewidth",1);
                clabel(C,h,[34 34.3 34.6],"color","k","labelspacing",40e3);
                if ~isempty(bsect)
                    patch([d/1e3; flip(d/1e3)],[bsect; flip(ssect)],"w");
                end
                if ii~=start
                    plot([d_init/1e3; flip(d_init/1e3)],[bsect_init; flip(ssect_init)],"-m");
                end
                P=patch([d(1)/1e3; d/1e3; d(end)/1e3],[-1500; Bsect; -1500 ],"w");
                P.FaceColor=[0.85 0.85 0.85];
                %colormap(gca,jet(10)); 
                %colormap(gca,othercolor("BuOrR_14"));
                Tmin = -1.4; Tmax = 1.3;
                %colormap(jet);
                dCM = ceil(abs(Tmin)*10);
                CM1 = othercolor("BuOrR_14",2*dCM); CM1 = CM1(1:dCM,:);
                dCM = ceil(3/4*abs(Tmax)*10);
                CM2 = othercolor("BuOrR_14",2*dCM); CM2 = CM2(dCM:end,:);
                dCM = ceil(abs(Tmax)*10)-dCM;
                CM3 = [linspace(CM2(end,1),152/255,dCM+1)' linspace(CM2(end,2),68/255,dCM+1)' linspace(CM2(end,3),158/255,dCM+1)'];
                CM = [CM1;CM2;CM3(2:end,:)]; colormap(gca,CM);
                caxis([Tmin Tmax]);
    
                ylim([-1.5e3 700]);
                xlim([0 max(d)/1e3-1]);
                ylabel("z [m]");
                xlabel("Distance [km]");
                grid on;
                box on;
    
                title("Temperature and Salinity");
    
                cb=colorbar(gca);
                cb.Location="southoutside";
                cb.Position=[0.11 0.6 0.22 0.012];
                cb.AxisLocation="out";
                set(cb,"Ticks",[-1.4:0.4:1.3]);
                set(cb,"TickLabels",{"-1.4","-1","-0.6","-0.2","0.2","0.8","1.2"});
                cb.XColor="k";
                cb.YColor="k";
                cb.TickLength=0.03;
                cb.FontSize=14;
                cb.Label.String = ["$\rm{T [}^\degree\rm{C]}$"];
                cb.Label.Interpreter = "latex";
    
                subplot("position",[0.37 0.75 0.3 0.2]); hold on;
    
                yyaxis left
                yticklabels({});
    
                yyaxis right
                plot(d_init/1e3,usect_init,"-m");
                plot(d/1e3,usect,"-k");
                xlim([0 max(d)/1e3-1]); xticklabels({});
                ylim([0 5e3]); ylabel("Ice speed [m/yr]");
                grid on; box on;
    
                text(0.5,1.15,["Time: ",titletime{kk}],"horizontalalignment","center","units","normalized","fontsize",18);
    
                subplot("position",[0.37 0.1 0.3 0.595]); hold on;
        %         pcolor(Dgrid,Zgrid,UVEL_parsect"); shading flat;
        %         UVEL_parsect(UVEL_parsect==0)=NaN;
        %         [C,h]=contour(Dgrid,Zgrid,UVEL_parsect",[0 0],"color",[1 1 1],"linewidth",1);
        %         if ~isempty(Ualonsect)
        %             patch([d(IsectinUa==1)/1e3 flipdim(d(IsectinUa==1)/1e3,2)],[bsect flipdim(ssect,2)],"w");
        %         end
        %         if ii~=2
        %             plot([d_init/1e3 flipdim(d_init/1e3,2)],[bsect_init flipdim(ssect_init,2)],"-m");
        %         end
        %         P=patch([d(1)/1e3 d/1e3 d(end)/1e3],[-1500 Bsect -1500 ],"w");
        %         P.FaceColor=[0.85 0.85 0.85];
        %         %colormap(gca,othercolor("RdBu7")); 
        %         colormap(gca,parula(10));
        %         caxis([-0.3 0.3]);
        % 
        %         ylim([-1.5e3 700]);
        %         xlim([0 max(d)/1e3-1]);
        %         ylabel(""); set(gca,"yticklabels",{});
        %         xlabel("Distance [km]");
        %         grid on;
        %         box on;    
        % 
        %         title("Speed parallel");
        % 
        %         cb=colorbar(gca);
        %         cb.Location="southoutside";
        %         cb.Position=[0.42 0.8 0.22 0.012];
        %         cb.AxisLocation="out";
        %         set(cb,"Ticks",[-0.3:0.1:0.3]);
        %         set(cb,"TickLabels",{"-0.3","-0.2","-0.1","0","0.1","0.2","0.3"});
        %         cb.XColor="k";
        %         cb.YColor="k";
        %         cb.TickLength=0.03;
        %         cb.FontSize=14;
        %         cb.Label.String = "Speed [m/s]";
    
                pcolor(Dgrid,Zgrid,VEL_perpsect); shading flat;
                VEL_perpsect(VEL_perpsect==0)=NaN;
                %[C,h]=contour(Dgrid,Zgrid,VEL_perpsect,[0 0],"color",[1 1 1],"linewidth",1);
                if ~isempty(bsect)
                    patch([d/1e3; flip(d/1e3)],[bsect; flip(ssect)],"w");
                end
                if ii~=start
                    plot([d_init/1e3; flip(d_init/1e3)],[bsect_init; flip(ssect_init)],"-m");
                end
                P=patch([d(1)/1e3; d/1e3; d(end)/1e3],[-1500; Bsect; -1500 ],"w");
                P.FaceColor=[0.85 0.85 0.85];
                %colormap(gca,othercolor("RdBu7")); 
                colormap(gca,parula(10));
                caxis([-0.1 0.1]);
    
                ylim([-1.5e3 700]);
                xlim([0 max(d)/1e3-1]);
                ylabel(""); set(gca,"yticklabels",{});
                xlabel("Distance [km]");
                grid on;
                box on;
    
                cb=colorbar(gca);
                cb.Location="southoutside";
                cb.Position=[0.42 0.6 0.22 0.012];
                cb.AxisLocation="out";
                set(cb,"Ticks",[-0.1:0.05:0.1]);
                set(cb,"TickLabels",{"-0.1","-0.05","0","0.05","0.1"});
                cb.XColor="k";
                cb.YColor="k";
                cb.TickLength=0.03;
                cb.FontSize=14;
                cb.Label.String = "Speed [m/s]";
    
                title("Speed perpendicular");
    
                subplot("position",[0.68 0.1 0.3 0.595]); hold on;
    
                %HeatFluxplot = HeatFlux_sect;%(:,:,kk);
        %        HeatFluxplot(HeatFluxplot==0)=999;
                %pcolor(Dgrid,Zgrid,HeatFluxplot); shading flat;

                TminTfplot = DelTheta_sect_in;
                pcolor(Dgrid,Zgrid,TminTfplot); shading flat;
    
                %[C,h]=contour(Dgrid,Zgrid,HeatFluxplot,[0 0],"color",[0.75 0.75 0.75],"linewidth",1);
                if ~isempty(bsect)
                    patch([d/1e3; flip(d/1e3)],[bsect; flip(ssect)],"w");
                end
                if ii~=start
                    plot([d_init/1e3; flip(d_init/1e3)],[bsect_init; flip(ssect_init)],"-m");
                end
                P=patch([d(1)/1e3; d/1e3; d(end)/1e3],[-1500; Bsect; -1500 ],"w");
                P.FaceColor=[0.85 0.85 0.85];
                colormap(gca,flipdim(othercolor("RdBu7"),1)); 
                %colormap(gca,parula(10));
                %caxis([-1e-2 1e-2]); %% FOR HEAT FLUX
                caxis([0 4]);
    
                ylim([-1.5e3 700]);
                xlim([0 max(d)/1e3-1]);
                ylabel(""); set(gca,"yticklabels",{});
                xlabel("Distance [km]");
                grid on;
                box on;
    
                cb=colorbar(gca);
                cb.Location="southoutside";
                cb.Position=[0.73 0.6 0.22 0.012];
                cb.AxisLocation="out";
                set(cb,"Ticks",[-0.01:0.005:0.01]);
                set(cb,"TickLabels",{"-0.01","","0","","0.01"});
                cb.XColor="k";
                cb.YColor="k";
                cb.TickLength=0.03;
                cb.FontSize=14;
                cb.Label.String = "Heat flux [TW]";
    
                title("Heat flux");
    
                pos = get(Hfig,"Position");
                set(Hfig,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
                fname = [folder,"/TSUQMsection_PIG_",num2str(runID),"_Time_",num2str(ii-1)];
                %print(Hfig,fname,"-dpng","-r400");
    
                hold off;
    
                currFrame = getframe(gcf);
                writeVideo(vidObj,currFrame);
    
                clf;
            end

        kk=kk+1;
    end

end

start = numel(nameFolds)+1;

if makemovie
    close(vidObj);
end

save("HeatVolumeTransport_"+figurefiletag+"_"+runID,"MITTime","DI","section","integral2D","dn","-v7.3");

return

%% some more plotting

[depth_m,MITTime_m] = ndgrid(-z,MITTime);

H=fig("units","inches","width",80*12/72.27,"fontsize",16,"font","Helvetica","defaultAxesColorOrder",[0 0 0; 0 0 0]);
subplot("position",[0.1 0.1 0.75 0.85]); hold on;

pcolor(MITTime_m,depth_m,HeatFlux_integralalongsect); shading flat;
g=plot(MITTime,bmean,"-w","linewidth",2);
datetick("x","mm/yyyy"); ylabel("Depth [m]");
%title("Heat transport");
legend(g,"Mean ice draft along section");
cb=colorbar; cb.Label.String="Heat transport through section [TW]";

ylim([-1200 -150]);

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = folder+"/Hovmoller_"+figurefiletag+"_"+runID;
print(H,fname,"-dpng","-r400");


H=fig("units","inches","width",80*12/72.27,"fontsize",16,"font","Helvetica","defaultAxesColorOrder",[0 0 0; 1 0 1]);
subplot("position",[0.1 0.1 0.75 0.85]); hold on;

yyaxis left 

plot(MITTime,[HeatFlux(:).sect_tot],"-k"); 
Mean1 = movmean([HeatFlux(:).sect_tot],12);
g(1)=plot(MITTime(7:end-7),Mean1(7:end-7),"-k","linewidth",2);
plot(MITTime,[HeatFlux(:).sect_in_tot],"-k");
Mean2 = movmean([HeatFlux(:).sect_in_tot],12);
g(2)=plot(MITTime(7:end-7),Mean2(7:end-7),"--k","linewidth",2);
plot(MITTime,[HeatFlux(:).sect_belowz],"-b");
plot(MITTime,[HeatFlux(:).sect_abovez],"-m");
%g(3)=plot(MITTime,HeatFlux_bt_sect,"-","color","r");
%g(4)=plot(MITTime,HeatFlux_bc_sect,"--","color","r");
%ylim([-5 6]); yticks([-5:6]); yticklabels({"","","","","-1","0","1","2","3","4","5","6"});
ylabel("Heat flux [TW]");
ylimits = ylim;

yyaxis right
plot(MITTime,Melt,"-","color",[0.7 0 1]);
Mean3 = movmean(Melt,12);
g(5)=plot(MITTime(7:end-7),Mean3(7:end-7),"-k","linewidth",2,"color",[1 0 1]);
% ylim([70 180]); yticks([70:10:180]); yticklabels({"70","80","90","100","110","","","","","","",""});
ylabel("Melt flux [TW]");
ylim(ylimits);

datetick("x","mm/yyyy"); 
grid on; box on;
%legend(g(:),{"Net heat transport","Heat inflow","Barotropic heat transport",...
%    "Baroclinic heat transport","Ice shelf melt flux"},"location","northwest");

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = [folder,"/HeatTransport_",figurefiletag,"_",runID];
print(H,fname,"-dpng","-r400");





    
