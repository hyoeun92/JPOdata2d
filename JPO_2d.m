steps = [1];

org=organizer('repository','./Models_Nakayama','prefix', 'PineIsland.','steps',steps);
md.miscellaneous.name='PineIsland';






%% {{{ Melt_param 1:
if perform(org,'Melt_param')

   md=loadmodel(org,'Parameterize');
   
	%Define constants
   g     = 9.81;
   E0    = 3.6e-2;
   Cd    = 2.5e-3;
   lam1  = -5.73e-2;
   lam2  = 8.32e-2;
   lam3  = 7.61e-4;
   Bs    = 7.86e-4;
   Bt    = 3.87e-5;
   CdTS0 = 6e-4;
   Ci    = 2.009*10^3;
   Cw    = 3.974*10^3;
   Ti    = -25;
   Li    = 3.35*10^5;
   Tm =  (Ci/Cw)*Ti - Li/Cw;
   
	 
   %Temporal arrays for the ambient temp, sal
   Toc = 0.47.*ones(md.mesh.numberofelements,1);
   Soc = 34.73.*ones(md.mesh.numberofelements,1);
   

   %Temporary
	md.mask.ice_levelset = -1.*ones(md.mesh.numberofvertices,1);
   pos = find(md.geometry.base==0);
   md.mask.ice_levelset(pos) = 1;
   icemask = mean(md.mask.ice_levelset(md.mesh.elements),2);
	oceanmask = mean(md.mask.ocean_levelset(md.mesh.elements),2);
   base = mean(md.geometry.base(md.mesh.elements),2);

	%load corrected base data to construct ambient stratified temperature, salinity
	load ./Models_JPO/base_corrected
	%Stratification
   strati_flag = true;
   
	if strati_flag == false

      Tdata=['T.0000004320.data'];
      Sdata=['S.0000004320.data'];
	   NOQ0=['SHICE_fwFlux.0000254880run23.data'];
		QQ0=['SHICE_fwFlux.0000254880run24.data'];
		TENQ0 = ['SHICE_fwFlux.0000254880run26.data'];
		nx = 460; ny=600; nz=130;
      Tdata1=readbin(Tdata,[460 600 130]);
      Sdata1=readbin(Sdata,[460 600 130]);
		noq01= readbin(NOQ0,[460 600]);
		qq01 = readbin(QQ0,[460 600]);
		tenq01=readbin(TENQ0,[460 600]);

	   xgOrigin = -102.838461072;
      ygOrigin = -75.4340;
      delX=0.007692308;
      delY=0.0020;
      delZ=10.0;
      xgrid = zeros(nx,1);ygrid = zeros(ny,1);zgrid=zeros(nz,1);
      xgrid(1,1) = xgOrigin;
      for i=2:nx
          xgrid(i) = xgrid(i-1) + delX;
      end
      ygrid(1,1) = ygOrigin;
      for i=2:ny
          ygrid(i) = ygrid(i-1) + delY;
      end
	   zgrid(1) = 0; %surface
      for i=2:nz
          zgrid(i) = zgrid(i-1) - delZ;
      end
      depth = zgrid;

	   lat = zeros(nx,ny); lon = zeros(nx,ny);
      for i=1:ny
          lon(:,i) = xgrid;
      end
      for i=1:nx
          lat(i,:) = ygrid;
      end
      [x1 y1]=ll2xy(lat,lon,-1,0,71);

		%interpolate melt from oceandata to ISSM grid
		noq02 = zeros(md.mesh.numberofvertices,1);
		tenq02 = zeros(md.mesh.numberofvertices,1);
      noq02 = griddata(double(x1(:)),double(y1(:)),double(noq01(:)),md.mesh.x,md.mesh.y,'nearest');
      qq02 = griddata(double(x1(:)),double(y1(:)),double(qq01(:)),md.mesh.x,md.mesh.y,'nearest');
		tenq02 = griddata(double(x1(:)),double(y1(:)),double(tenq01(:)),md.mesh.x,md.mesh.y,'nearest');

		%meltrates from oceanmodeldata is negative, so convert it to be positive
		tenq02 = -tenq02 ;  noq02 = -noq02; qq02 = -qq02;
		tenq02 = (365*24*3600*tenq02)/(md.materials.rho_ice);
      noq02 = (365*24*3600*noq02)/(md.materials.rho_ice); 
      qq02 = (365*24*3600*qq02)/(md.materials.rho_ice);
 
      %Choose layers to interpolate
	   first = 30; endd = 89; %
      depth_chosen = depth(first:endd);
	   Tdata2 = zeros(md.mesh.numberofvertices,length(depth_chosen)); Sdata2 = zeros(md.mesh.numberofvertices,length(depth_chosen));
    
	   for k=first:endd
          kk = k - (first-1);
          Tdata2_tmp = Tdata1(:,:,k);
          Tdata2(:,kk) = griddata(double(x1(:)),double(y1(:)),double(Tdata2_tmp(:)),md.mesh.x,md.mesh.y,'nearest');
          Sdata2_tmp = Sdata1(:,:,k);
          Sdata2(:,kk) = griddata(double(x1(:)),double(y1(:)),double(Sdata2_tmp(:)),md.mesh.x,md.mesh.y,'nearest');
      end
      temp_input = Tdata2; sal_input = Sdata2;
   
	   for k=1:length(depth_chosen)
          pos00 = find(temp_input(:,k) ==0);
          temp_input(pos00,k) = mean(nonzeros(temp_input(:,k)));
          pos01 = find(sal_input(:,k) <=32); pos02 = find(sal_input(:,k)>32);
          sal_input(pos01,k) = mean(nonzeros(sal_input(pos02,k)));
      end
      %Interpolate temperature and salinity fields
	   temp_atbase = zeros(size(md.geometry.base)); sal_atbase = zeros(size(md.geometry.base));
      temp_depths = depth_chosen;
  
	   %binary search
      base_plus = -md.geometry.base; temp_depths_plus = -temp_depths; %sal_depths_plus = -sal_depths;

      for i = 1:md.mesh.numberofvertices
          if base_plus(i)<0
             temp_atbase(i) =0;
             sal_atbase(i) = 0;
          else
             offset = binaryy(base_plus(i),temp_depths_plus,length(temp_depths));
             if offset ==-1
                temp_atbase(i) = temp_input(i,1);
                sal_atbase(i)  = sal_input(i,1);
             elseif offset == length(temp_depths)
                temp_atbase(i) = temp_input(i,end);
                sal_atbase(i)  = sal_input(i,end);
             else
                deltaz = temp_depths_plus(offset+1) - temp_depths_plus(offset);
                alpha2 = ( base_plus(i) - temp_depths_plus(offset) )/deltaz;
                alpha1 = 1.0 - alpha2;
                temp_atbase(i) =  alpha1*temp_input(i,offset) + alpha2*temp_input(i,offset+1);
                sal_atbase(i)  =  alpha1*sal_input(i,offset) + alpha2*sal_input(i,offset+1);
             end
          end
      end %binary search loop
  
       save ./Models_JPO/stratified_amibient temp_atbase sal_atbase Tdata2 Sdata2 temp_input sal_input	
       save ./Models_JPO/oceanmodelmelt noq02 qq02 tenq02
	 else %stratiflag
       load ./Models_JPO/stratified_ambient.mat
       load ./Models_JPO/oceanmodelmelt.mat 
   end
	
   %% Calculate grounding line depth (Pelle et al.(2023))
   zgl_cal = true;
	if zgl_cal == false;
	
     %Get reldistgl
     filename = tempname();
     contours = isoline(md,md.mask.ocean_levelset,'value',0) ;
     expwrite(contours,'yam.exp');
     distgl = exp2levelsetfunction(md,'yam.exp');
     distgl_temp = distgl;
     distgl = abs(mean(distgl(md.mesh.elements),2));
     filename = tempname();
     contours=isoline(md,md.mask.ice_levelset,'value',0);
     expwrite(contours,'yam2.exp');
     distif = exp2levelsetfunction(md,'yam2.exp').*-1;
     distif_temp = distif;
     distif = abs(mean(distif(md.mesh.elements),2));
     reldist_gl = distgl./(distgl + distif);

	  %Interpolate fields and get base of ice shelf
     Xq = min(md.mesh.x):300:max(md.mesh.x);
     Yq = max(md.mesh.y):-300:min(md.mesh.y);
     [X Y]=meshgrid(Xq,Yq);
     B = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.geometry.base,Xq,Yq,NaN);
     maskg = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,Xq,Yq,NaN);
     maski = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,Xq,Yq,NaN);
     reldist = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,reldist_gl,Xq,Yq,NaN); %relative distance interpolated onto XQYqGRID
     B(find(maskg>0 | maski>0)) = NaN; %Only want ice shelf 
   
	  % Turn the head dem into a GRID object
     DEM    = GRIDobj(X,Y,B); % Creates a DEM object from the vector(coordinate) X and Y and the matrix dem, The elements of dem refer to the elevation of eah pixel.
     DEMf   = fillsinks(DEM); % removes topographi depressions.
     FD     = FLOWobj(DEMf); % FD contains a topological ordering of the DEMF nodes(FLOWobj creates a flow direction object)
     [U,V]  = flowvec(FD); %flow velocity vector from FLOWobj
   
	  %Get drainage basins
     D = drainagebasins(FD); %drainagebasins(FD) calculates the drainage basins from all drainage network outlets
   
	  %Compute minimum point in drainage basin
     zglgrid = B; % md.geometry.base that is interpolated onto Xq Yq grid
     %Loop over drainage basins
     for i=min(D.Z(:)):max(D.Z(:))
        %Get cells within drainage basin i, then get reldistgl and minimum base depth
        pos = find(D.Z(:)==i);         % Get cells within drainage basin i
        reldistbasin = reldist(pos);   % reldistgl  reldist is relative distance from the GL that is interpolated onto XqYq GRID
        [zmin loc] = min(DEM.Z(pos));  % minimum basin depth
        reldistpt(i+1) = reldistbasin(loc); %reldistgl at point of minimum base depth
          if abs(reldistpt(i+1))<0.3        %check that basin drains over grounding line and not ice front
             zglgrid(pos) = zmin;           %set zgl of drainage basin i equal to min base depth
          end
        clear pos;
     end

	  %Interpolate back onto mesh
     D.Z=cast(D.Z,'double');
     DB = InterpFromGridToMesh(Xq',fliplr(Yq)',flipud(D.Z),md.mesh.x,md.mesh.y,NaN);
     pos=find(DB~=0);
     pos2=find(DB==0);
     DB(pos2) = round(griddata(md.mesh.x(pos),md.mesh.y(pos),DB(pos),md.mesh.x(pos2),md.mesh.y(pos2),'nearest'));
     zgl = InterpFromGridToMesh(Xq',fliplr(Yq)',flipud(zglgrid),md.mesh.x,md.mesh.y,NaN); %%%%%%
     zgl_vertices = zgl;
	  pos = find(md.geometry.base==0);
     zgl_vertices(pos) = 0;
     zgl = mean(zgl_vertices(md.mesh.elements),2);
     clear pos;

	  save ./Models_JPO/grounding_line_depth zgl_vertices zgl 
     else
     
     load ./Models_JPO/grounding_line_depth
     
     end
	  
	  %Maximal basal slope from corrected base data
     [sx,sy,s] = slope(md,base_corr);
     alpha = atan(sqrt(sx.*sx + sy.*sy));
     K = (1/32)*ones(4);  % 1/16 for nonq0 1/32 for Q0, 1/64 for 10Q0 
     alpha = conv2(alpha,K,'same') ; %smoothing
     
	  % Load discharge data
     % Discharge data is derived from running a thermal model. Details can be found in Nakayama et al.(2021)
	  % The following locations are where the data is defined. 
     dis_den = zeros(md.mesh.numberofvertices,1);
	  
	  %R1 in Nakayama et al.(2021) 
	  locR1 = 14768;
	  locR2 = 21392;
	  dis_den(locR1) =  200*20*0.011; 
     md.mask.ocean_levelset(locR1)=0;
	  %R2 in Nakayama et al.(2021) 
     dis_den(locR2) =  200*20*0.0036; 
	  md.mask.ocean_levelset(locR2)=0;

	  
     dis_den_temp = dis_den;

	  %Get discharge along primary grounding line
     filename     = 'yam33.exp' ;%tempname();
     contours     = isoline(md,md.mask.ocean_levelset,'value',0);
     allVals      = [contours.nods];
     [maxval pos] = max(allVals); %for primary grounding line
     contours     = contours(pos);% contours for primary grounding line
     expwrite(contours,'yam33.exp');% create expfile for the primary grounding line contours
     [elements,x,y,z1,horz,disgl]=SectionValues(md,dis_den,'yam33.exp',[50 1]); %get the value of  discharge data following the PGL contours
     
	  compareNAKA_flag = true;
     
	  if compareNAKA_flag ==false % We don't use channnel width algorithm here, as channel width of 200 m is prescribed in Nakayama et al. (2021)

	    %Find channel widths for channels discharging >2 m3/s
       pos2 = find((disgl)>2);
	    %check to make sure we have a channel width
       disflg = false;
       %If we do not have discharge draining across the GL, use initial glfw conditions
       if  length(pos2)<2
           disflg = true;
           disp('Value of FW is below 2m^{3}s{-1}')
      
       else %If we have discharge across the GL, compute channel width
           difference = diff(pos2); %size is dim(pos-1,1), ex) for pos2 = [3 4 5], diff(pos2) =[1 1]
           b = find([difference' inf]>1); %find position difference larger than 1 from the matrix with size is dim(pos2,1) that contains inf
           c = diff([0 b]);
           d = cumsum(c); %cumulative sum
           groups = [];
             for i=1:length(d)
                if i==1
                %Compute distance from first to last GL discharge point
                x_diff(i) = x(pos2(d(i)))-x(pos2(1));
                y_diff(i) = y(pos2(d(i)))-y(pos2(1));
                groups(i) = sqrt(x_diff(i)^2+y_diff(i)^2);
                else
                %Compute distance from first to last GL discharge point
                x_diff(i) = x(pos2(d(i)))-x(pos2(d(i-1)+1));
                y_diff(i) = y(pos2(d(i)))-y(pos2(d(i-1)+1));
                groups(i) = sqrt(x_diff(i)^2+y_diff(i)^2);
                end
             end
            groups(find(groups<50))=[];
	         if isempty(groups)
            disflg = true;
            disp('Channel width cannot be defined as it is smaller than 50m')
            %chan_wid = initdis_params{1};
            else
            chan_wid = mean(groups);
            end
       end %length(pos2) if

	 else
		  disflg = false;
        chan_wid = 200; %for the comparsion of Nakayama et al.(2021)
	  end % if statements for compareNAKA_flag ends
	 	  
	%Apply channel width to discharge flux (scale from m3/s to m2/s)
   dis_den   = dis_den./chan_wid; %% size(md.mesh.numberofvertices)
	
	% maxdist is originally set to be 5000 in Pelle et al.(2023), but it is changed to be 30 
	% This is because of the calculation of the point-like discharge fields defined in Nakayama et al.(2021) 
   maxdist   = 30; 

   %Set discharge over floating ice to 0
   dis_den_gr = dis_den;  % Q_{0} = ? m^{2}s{-1}
   dis_den_gr(find(md.mask.ocean_levelset<0))=0;
   
   %Calculate distance from grounding line discharge location field
   filename = tempname();
   contours=isoline(md,dis_den_gr,'value',10^-9);  % originally, (2/chan_wid) creating contour where nonzeroQ0 larger than 2/chan_wid &  grounded elements
                                                   % For the comparison with Nakayama et al(2021) where the amount of fw was prescribed, don't have to use 2/chan_wid
   expwrite(contours,'yam3.exp');
   dist_disden = exp2levelsetfunction(md,'yam3.exp'); %yam3.exp are two contours surrounding two discharge outflow (point-source like) locations
   
	%Create mask where discharge will be applied (1 where discharge, 0 elsewhere)
   dismask = zeros(size(dist_disden));
	pos4 = find(dist_disden<=maxdist); % for pointwise search
   dismask(pos4) = 1;  %dismask is 1 over floating elements where radial distance is < maxdist (where initial Q0 will be applied)
   
	%Create contours around the mask so we can access elements within
   filename = 'yam4.exp';  %tempname();
   contours = isoline(md,dismask,'value', 0.0001);  %creating contour around the region where dismask = 1
   expwrite(contours,'yam4.exp');
   profiles=expread('yam4.exp');
   
	%Additional treatment for the point-like discharges defined in Nakayama et al. (2021)
   dismask2 = zeros(size(dist_disden));
	maxdist2 = 5000;
   pos44 = find(dist_disden<=maxdist2);
   dismask2(pos44) = 1;
	contours2 = isoline(md,dismask2,'value', 0.0001);  %creating contour around the region where dismask = 1
   expwrite(contours2,'yam44.exp');
   profiles2=expread('yam44.exp');
	profile_1 = profiles2(1) ; 
	profile_2 = profiles2(2) ; 
   ele_1   = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,profile_1,'element',0)); 
   ele_2   = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,profile_2,'element',0));   
   
	dis_den_elements = dis_den(md.mesh.elements);
   oceanmask_elements = md.mask.ocean_levelset(md.mesh.elements);
   icemask_elements = md.mask.ice_levelset(md.mesh.elements);
	temp_atbase_elements = temp_atbase(md.mesh.elements); 
	sal_atbase_elements = sal_atbase(md.mesh.elements);

	[ele row] = find(dis_den_elements~=0);
   max_dis_den_elements = max(dis_den_elements(ele,1:3),[],2);
	max_temp_elements    = max(temp_atbase_elements(ele,1:3),[],2); 
	max_sal_elements     = max(sal_atbase_elements(ele,1:3),[],2);

   for l=1:length(ele)
   dis_den_elements(ele(l),1:3) = max_dis_den_elements(l);
   oceanmask_elements(ele(l),1:3) = 0;
	icemask_elements(ele(l),1:3) = -1;
   temp_atbase_elements(ele(l),1:3) = max_temp_elements(l) ;
   sal_atbase_elements(ele(l),1:3) =  max_sal_elements(l)  ;
	end
   dis_den = mean(dis_den_elements,2);
   oceanmask = mean(oceanmask_elements,2);
   icemask = mean(icemask_elements,2);
	Toc = mean(temp_atbase_elements,2);
	Soc = mean(sal_atbase_elements,2) ;
   
	plotmodel(md,'data',base_corr,'figure',1001);
	posg_vertices = find(md.mask.ice_levelset ~= -1);
	plotmodel(md,'data',base_corr,'figure',1002);
 
	boundaryice1 = find(icemask>-1 & icemask<1 );
	icemask(boundaryice1) = 1;   
   
	boundaryice2 = find(oceanmask>-1 & oceanmask<0); 
   icemask(boundaryice2) = 1;

	boundaryice3 = find(oceanmask>0 & oceanmask<=1);
   icemask(boundaryice3) = 1;
   
	isice_ele_1 = find(oceanmask(ele_1)>-1 & oceanmask(ele_1)<=0 );
   icemask(ele_1(isice_ele_1)) = -1;
   isice_ele_1 = find(oceanmask(ele_1)>=0 & oceanmask(ele_1)<1 );
   icemask(ele_1(isice_ele_1)) = -1;

	isice_ele_2 = find(oceanmask(ele_2)>-1 & oceanmask(ele_2)<=0 );
   icemask(ele_2(isice_ele_2)) = -1;
   isice_ele_2 = find(oceanmask(ele_2)>=0 & oceanmask(ele_2)<1 );
   icemask(ele_2(isice_ele_2)) = -1;

	% In general, floating ice shelf can be searched by icemask<0 & oceanmask<0
	% but the PIG domain employed in Nakayama et al(2021) -- icemask ==-1 stands for floating ice shelf where melting is applied
   posg = find(icemask~=-1);
   alpha(posg) = 0;
   Toc(posg) = 0; Soc(posg) = 0;
  
	%alpha treatment
	alpha_tr0   = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'alpha_treat0.exp','element',0)); 
   alpha_tr3   = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'alpha_treat3.exp','element',0));
	K_add = (1/16)*ones(4); 
	alpha(alpha_tr0) = conv2(alpha(alpha_tr0),K_add,'same') ;
   alpha(alpha_tr3) = conv2(alpha(alpha_tr3),K_add,'same') ;
	

	% Define grounding line freezing temperature having size(numberofelements,1) before loop
   T_ini = lam1*0. + lam2 + lam3.*zgl;  % freezing temperature of meltwater plumes at the grounding line: T_f(0,zgl)
   Fr_gl = lam1.*Soc + lam2 + lam3.*zgl;  % T_f(S_amb,zgl)
   
	%Define variables for the parameterization
   Ycd = zeros(md.mesh.numberofelements,1); %coordinate on which the parameterization is defined 
	Tf_z = zeros(md.mesh.numberofelements,1); %T_f(S_amb,z)
	M0_app = zeros(md.mesh.numberofelements,1);
	source_correction = zeros(md.mesh.numberofelements,1); 
	balance_tf_field_temp_gl = zeros(md.mesh.numberofelements,1);
   balance_tf_field_temp = zeros(md.mesh.numberofelements,1); 
	balance_tf_field_temp_weighted=  zeros(md.mesh.numberofelements,1);   
	balance_tf_field = zeros(md.mesh.numberofelements,1);
   balance_tf_field_with_corr = zeros(md.mesh.numberofelements,1); 
   U_ini_field_temp = zeros(md.mesh.numberofelements,1); 
	U_ini_field = U_ini_field_temp;
   U_ini_balancefield_temp = zeros(md.mesh.numberofelements,1);
	U_ini_balancefield = zeros(md.mesh.numberofelements,1);
   balance_vel_field_temp1 = zeros(md.mesh.numberofelements,1);
   balance_vel_field_temp2 = zeros(md.mesh.numberofelements,1);
   balance_vel_field = zeros(md.mesh.numberofelements,1);	
	sinalp = zeros(md.mesh.numberofelements,1);  
	alpha_u = zeros(md.mesh.numberofelements,1); 
	del_rho_m = zeros(md.mesh.numberofelements,1);
   distfield_tf = zeros(md.mesh.numberofelements,1);  
	distfield_vel = zeros(md.mesh.numberofelements,1);
   VEL_profile = zeros(md.mesh.numberofelements,length(profiles));
	TEMP_profile = zeros(md.mesh.numberofelements,length(profiles)); 
	MELT_profile = zeros(md.mesh.numberofelements,length(profiles)); 

	meltrate_field = zeros(md.mesh.numberofelements,1);
   

	pos = find(alpha>=pi) ; alpha(pos) = alpha(pos) - 0.001; clear pos;
   pos = find(zgl>base); zgl(pos) = base(pos); clear pos;

   %% Set ambient temperature bound (Lazerom et al. (2018))
   posf = find(icemask==-1);
   pos2 = find(Toc(posf)<=lam1.*Soc(posf) + lam2);
   pos3 = intersect(posf(pos2),posf);
   Toc(pos3) =lam1.*Soc(pos3) + lam2 ;
   

	clear posf; clear pos2; clear pos3;

	pos_minalpha = find(alpha<10^-20 & icemask==-1);  
   alpha(pos_minalpha) = 10^-20;
   sinalp = sin(alpha);

	%velocity parameter
   alpha_u = 0.5585.*(sinalp).^(-0.6858); 
	alpha_u_gl = 0.333; % for initial field description of 2D

	Tf_z = Toc - (lam1.*Soc + lam2 + lam3.*base ); % T_{ambs} - T_f(S_amb,zb)
   M0_app = CdTS0./( (lam1.*Soc + lam2 + lam3.*base  )  - Tm ); %meltrate factor defined in Jenkins (2011)
	del_rho_m = (Bs.*(Soc)-Bt.*(Toc - Tm)).*g.*sinalp; % buyancy contrast defined using effective meltwater temperature (Beckmann et al)

   % 1) Loop over contours, extract maximum discharge in contour, and apply to all elements within contour
   % 2) computing 5L' is not necessary, instead, computing velocity and thermal forcing fields
   
  if disflg
      disp('No discharge across the GL')
  else %If discharge across GL
      for i=1: size(profiles,2)
	       %Compute max discharge in each region
          profile=profiles(i);
          isdis_node = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,profile,'node',2));  %find node inside the contour
          isdis_el   = find(ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,profile,'element',0)); %find element inside the contour	  
          if (isempty(isdis_el))
               disp('no discharge inside the contour')
               maxdis(i) = 0. ;  
               ystar(i) =  0. ;

           else
               %Find maximum discharge close to GL
               alldis       = dis_den(isdis_el); %get averaged discharge values inside the profile
               allmask      = oceanmask(isdis_el); %get oceanmask values inside the profile
               pos          = find(alldis~=0 & allmask>=0); %find elements where nonzero fw near the GL (original : allmask>=0) 
               posmin       = find(allmask(pos)<min(allmask(pos)+40));  
               if isempty(posmin)
                  maxdis(i) = 0.;
                  ystar(i) = 0.;
               else
                  [maxdis(i) loc] = max(alldis(pos(posmin))); %a maximal fw  value  & location
                  dist_dis_temp = zeros(md.mesh.numberofvertices,1);
                  dist_dis_temp(md.mesh.elements(isdis_el(pos(posmin(loc))),:)) = 1; % it has one nonzero value  at  Q0_{max} near the GL 
                  contours=isoline(md,dist_dis_temp,'value',1);  %creating point contour
                  expwrite(contours,'yam5.exp');
                  dist_dis_temp2 = exp2levelsetfunction(md,'yam5.exp'); % calculating distances from the Q0_max; the distances are all positive as it is from the point contour  
       
                     
						%calculate the grounding line quantities that is averaged over maxdist
                  slopeavg = mean(nonzeros(sin(alpha(isdis_el))));
                  T_diff_ini = mean(abs(Toc(isdis_el)-T_ini(isdis_el)),'omitnan');
                  S_diff_ini = mean(abs(Soc(isdis_el)),'omitnan');
                  del_rho_gl = (Bs.*S_diff_ini-Bt.*T_diff_ini).*g;  %same
                  del_rho_m_gl = (Bs.*( mean(Soc(isdis_el),'omitnan') )-Bt.*( mean(Toc(isdis_el),'omitnan') - Tm )).*g.*slopeavg;
                  Tf_gl = mean( abs( Toc(isdis_el) - Fr_gl(isdis_el) ) ,'omitnan' )  ; % ambient thermal forcing averaged over the gl source loc
                  M0_app_gl_tmp = mean(Fr_gl(isdis_el),'omitnan');
                  M0_app_gl = CdTS0./(M0_app_gl_tmp - Tm) ;
                  Tc_gl = Tf_gl./(1. - (M0_app_gl*Tm)/(E0*slopeavg) ); %balance thermal forcing defined at the gl source loc
						GT_gl = 1./(1. - (M0_app_gl*Tm)/(E0*slopeavg) );

                  %Define finite source
                  ystar(i) =  ( (maxdis(i)^2.)./((E0.^2.).*(slopeavg.^3.).*del_rho_gl ) ).^(1.0/3.0);  
                  zstar(i) =  ystar(i)*slopeavg ; 

						%thermal_forcing_parameters
                  alpha_t = 1.2182*(GT_gl^(-0.72))*(Tc_gl^(-0.0231));
                  c1 = 0.0065; c2 = -0.0022; c3 = 0.67; c4 =-0.6635; c5  =-0.0486; c6 =0.0651;
                  r1 = 0.53119; r2 = 0.2553; r3= 0.2581 ; r4 = -1.;
                  K1 =1. ; K2 = 1. ;
                  A1 = 2.1449; A2 = -0.0094;
                  b1 = -0.1003;

						Km(i) = ( (K1*maxdis(i))^(c1*(Tc_gl)^r1 + c2) )*( (K2*Tc_gl)^( (c3*(maxdis(i)^r2)) + (c4*(maxdis(i)^r3)) + (c5*(Tc_gl^r4)) + c6   ) );
                  abs_zgl = mean(abs(zgl(isdis_el)),'omitnan');
                  final_coeffi(i) = (  A1*(abs_zgl^b1)*(1 + A2*Tf_gl ) )*Km(i);
 
						%Source correction parameter
                  corrected_depth(i) =  -(abs_zgl + zstar(i)); %Depth should be negative
						DTFDP = -0.0743; %dT_{f}/dP for fw, degree/MPa, Millero (1978)
                  T_melting_corr(i) = DTFDP*(1028*9.81*abs(corrected_depth(i))*(10^-6)-0.1);
						Tf_source(i) = T_melting_corr(i) - (lam2 + lam3*(corrected_depth(i))); % Thermal forcing defined at the virtual source location

                  Ycd = mean(dist_dis_temp2(md.mesh.elements),2); %coordinate x on which the parameterization is based on
                  posf = find(icemask==-1);  %find floating ice ;
               
						%Calculate initial velocity fields & initial balance velocity fields with decaying exponent 1/3
						U_ini_field_temp(posf) = ((maxdis(i)*del_rho_gl*slopeavg)/(Cd + E0*slopeavg))^(1./3.) ;  
                  U_ini_field(posf)  = U_ini_field_temp(posf) -  U_ini_field_temp(posf).*( 1. - (ystar(i)./(ystar(i) +  Ycd(posf))).^(alpha_u_gl )  ) ; %Eq. 50
                  U_ini_balancefield_temp(posf) = ((maxdis(i).*del_rho_gl.*sinalp(posf))./(Cd + E0*sinalp(posf))).^(1./3.) ; %u_{c0,\theta_{local}} in Eq.(41) 
                  U_ini_balancefield(posf)  = U_ini_balancefield_temp(posf) -  U_ini_balancefield_temp(posf).*( 1. - (ystar(i)./(ystar(i) +  Ycd(posf))).^(alpha_u_gl )  ) ; 

                  %Calculate thermal forcing fields
                  distfield_tf(posf) = (M0_app_gl)*(final_coeffi(i))*(1. - ( ystar(i)./(Ycd(posf) + ystar(i)) ).^(alpha_t));
                  balance_tf_field_temp(posf) = (Tf_z(posf))./(1 - (M0_app(posf).*Tm)./(E0*sinalp(posf))  ); % balance thermal forcing defined by local slope angle, sinalp
                  mean_angle = mean(sin(alpha),'omitnan');  % basal slope angle representative of the shelf cavity of interest
                  gl_angle  = slopeavg; % averaged-gl slope angle 
                  balance_tf_field_temp_gl(posf) = (Tf_z(posf))./(1 - (M0_app(posf).*Tm)./(E0*gl_angle )  ); % balance thermal forcing defined by the grounding line slope angle, slopeavg
                  epsilon_mean = (E0*mean_angle)/(CdTS0 + E0*mean_angle); %thermal forcing efficiency defined by mean canvity angle(Jenkins (2018))
                  epsilon_gl  = (E0*gl_angle )/(CdTS0 + E0*gl_angle ); %thermal forcing efficiency defined by the grounding line slope angle(Jenkins (2018))
                  epsilon_coeffi=2.;
                  weight1 = epsilon_gl^epsilon_coeffi ;  
						weight2 = epsilon_mean^epsilon_coeffi;
                  
                  balance_tf_field_temp_weighted(posf)  = ( weight1.*balance_tf_field_temp(posf) + weight2.*balance_tf_field_temp_gl(posf) )./(weight1 + weight2); %weighted balance thermal forcing fields 
                  balance_tf_field(posf)  = (M0_app_gl)*(final_coeffi(i)).*(balance_tf_field_temp_weighted(posf)).*(1. - ( ystar(i)./(Ycd(posf) + ystar(i)) ).^(alpha_t));
                  source_correction(posf) = (M0_app_gl)*(final_coeffi(i))*(Tf_source(i))*(1. - ( ystar(i)./(Ycd(posf) + ystar(i)) ).^(alpha_t));
                  balance_tf_field_with_corr(posf) = balance_tf_field(posf) + source_correction(posf); 

                  %Calculate velocity fields
                  distfield_vel(posf) = (1. - ( ystar(i)./(Ycd(posf) + ystar(i)) ).^(alpha_u(posf)));
                  balance_vel_field_temp1(posf) = (  M0_app(posf).*(2.0/3.0).*(del_rho_m(posf)./(Cd + E0.*sinalp(posf)) ).*(Ycd(posf))).*(Tc_gl.*((ystar(i)./(ystar(i) +  Ycd(posf))).^( alpha_u_gl ) ));
						balance_vel_field_temp2(posf) = sqrt( U_ini_balancefield(posf).^2. + balance_vel_field_temp1(posf)); %sqrt of Eq.40
                  balance_vel_field(posf) = (U_ini_field(posf)) + (balance_vel_field_temp2(posf) - U_ini_field(posf)).*(distfield_vel(posf)); % Eq.49

						% save for ith profile
                  TEMP_profile(posf,i) = balance_tf_field_with_corr(posf);
						VEL_profile(posf,i) = balance_vel_field(posf);


				end %inner if else
			 end %outer if else
        end %for profiies loop

		    
	  end %if else disflg 
     
	      
  
  posf = find(icemask==-1);  %for floating elements
  
  % MELT calculation based on U\Delta T, note that each fields driven by each fw are superimposed
  MELT_profile(posf,1) = TEMP_profile(posf,1).*VEL_profile(posf,1);
  MELT_profile(posf,2) = TEMP_profile(posf,2).*VEL_profile(posf,2);
  meltrate_field(posf) = 0.5.*(MELT_profile(posf,1) + MELT_profile(posf,2)) ; 
  
  meltrate_field_year = 365*24*3600.*meltrate_field;

  plotmodel(md,'data',meltrate_field_year,'mask', icemask==-1,'caxis',[0 170],'figure',101,'colorbartitle','Parameterized melt rates (m yr^{-1})','xticklabel', [],'yticklabel',[],...
  'axis','off','colorbarcornerposition',  'south', 'colormap',slanCM(12)); 
  
  % For ocean model melt rates
  noq02  = mean(noq02(md.mesh.elements),2); 
  qq02   = mean(qq02(md.mesh.elements),2);
  tenq02 = mean(tenq02(md.mesh.elements),2);
  
  plotmodel(md,'data',qq02,'mask', icemask==-1,'caxis',[0 170],'figure',103,'colorbartitle','Melt rates from the ocean model (m yr^{-1})','xticklabel', [],'yticklabel',[],'axis','off',...
  'colorbarcornerposition', 'south','colormap',slanCM(12));
  plotmodel(md,'data',noq02,'mask', icemask==-1,'caxis',[0 170],'figure',104,'colorbartitle','Melt rates from the ocean model (m yr^{-1})','xticklabel', [],'yticklabel',[],'axis','off',...
	   'colorbarcornerposition', 'south','colormap',slanCM(12)); % 37 13 9 12 
  plotmodel(md,'data',tenq02,'mask', icemask==-1,'caxis',[0 250],'figure',105,'colorbartitle','Melt rate (m/yr)','xticklabel', [],'yticklabel',[],'axis','off',...
      'colorbarcornerposition', 'south','colormap',slanCM(12)); % 37 13 9 12

end   %org
%% }}}	
