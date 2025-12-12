% Viscous hydro-thermomechanical-chemical (HTMC) code
% Solving momentum, mass and energy conservation eqs.
% for coupled liquid-solid system (no hydrous melting)
% in primitive variable formulation
% with variable viscosity and thermal conductivity
% using FD with staggered grid and marker in cell techniques

% BEFORE RUNNING:
% SELECT VENUS COMPOSITION
% SELECT FLOW LAW
% SET GEOMETRY AND INITIAL CONDITIONS

% Clearing memory and figures
clf
clear all;% clf

% --- VIDEO OPTIONS ---
%videoFile = 'plagiocase_video.mp4'; % Name of the video file
%v = VideoWriter(videoFile, 'MPEG-4'); % Use 'MPEG-4' for .mp4 file
%v.FrameRate = 5;  % For example, set the FPS to 30
%open(v); % Open the video for writing

% --- DEBUG OPTIONS ---
%dbstop if error          % stop on any error
%dbstop if naninf         % stop when NaN or Inf appears
% treat divide-by-zero as an error (instead of silent Inf)
%warning('error', 'MATLAB:divideByZero');

% Database | Hydration mol frac | + Hydration anomaly | + temp anomaly
nname ='VENERA-14_0.09_0.03_300_';  % mat filename
timestep=0;
% Load mat file
fdata=fopen('file.txt','rt');
timestep=fscanf(fdata,'%d',1);
fclose(fdata);
if(timestep>0)
namemat    =  [nname,num2str(timestep)];
load(namemat);
timestep=timestep+1;
else
% Start new run
output_dir = fullfile(pwd, nname);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Select composition of Venusian basalt from different missions (1-6 see below)
Venus_composition=5;
% 1: Kargel et al. (1993), VEGA
% 2: Surkov et al. (1986), VEGA-2
% 3: Semprich et al. (2024), Venusian Basalt
% 4: Filiberto et al. (2014); Fegley (2003), Venera-13
% 5: Filiberto et al. (2014); Fegley (2003), Venera-14
% 6: McDonough & Sun (1995), BSE
composition_names={'VEGA', 'VEGA-2', 'VBasalt','Venera-13', 'Venera-14', 'BSE'};

% Load all data from the database
namemat = sprintf('database_%s.mat', composition_names{Venus_composition});
filepath = fullfile('databases', namemat);
load(filepath);

% Define Numerical model
xsize = 40000;         % Horizontal model size, m
ysize = 25000;         % Vertical model size, m
Nx    = 81;            % Horizontal grid resolution
Ny    = 51;            % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

% Define Gravity
gx=0; % Horizontal gravity acceleration, m/s^2
gy=8.87; % Vertical gravity acceleration, m/s^2 (Venus)

% Coordinates of different nodal points
% Basic nodes
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
% Vx-Nodes
xvx=0:dx:xsize+dy; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
% Vy-nodes
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
% P-Nodes
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% Nodal arrays
% Basic nodes
ETA=zeros(Ny,Nx); % Viscosity, Pa*s
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
TYP=zeros(Ny,Nx); % Marker type     
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Total density, kg/m^3
RHOFX=zeros(Ny1,Nx1); % Fluid density, kg/m^3
KX=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIX=zeros(Ny1,Nx1); % Porosity
vx=zeros(Ny1,Nx1); % Solid vx-velocity m/s
vxf=zeros(Ny1,Nx1); % Fluid vx-velocity m/s
RX=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qxD=zeros(Ny1,Nx1); % qx-Darcy flux m/s
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Total density, kg/m^3
RHOFY=zeros(Ny1,Nx1); % Fluid density, kg/m^3
KY=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIY=zeros(Ny1,Nx1); % Porosity
vy=zeros(Ny1,Nx1); % Solid vy-velocity m/s
vyf=zeros(Ny1,Nx1); % Fluid vy-velocity m/s
RY=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qyD=zeros(Ny1,Nx1); % qy-Darcy flux m/s
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
RHOWATER=zeros(Ny1,Nx1); % Water density, kg/m^3
RHOSOLID=zeros(Ny1,Nx1); % Solid density, kg/m^3
RHOCP=zeros(Ny1,Nx1); % Volumetric heat capacity, J/m^3/K
ALPHA=zeros(Ny1,Nx1); % Adiabatic heating term
ALPHASOLID=zeros(Ny1,Nx1); % Thermal expansion solid, 1/K
ALPHAFLUID=zeros(Ny1,Nx1); % Thermal expansion fluid, 1/K
BETTASOLID=zeros(Ny1,Nx1); % Compressibility solid, 1/Pa
BETTAFLUID=zeros(Ny1,Nx1); % Compressibility  fluid, 1/Pa
HR=zeros(Ny1,Nx1); % Radioactive heating, W/m^3
HA=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMAxx, 1/s
tk1=zeros(Ny1,Nx1); % Old temperature solid, K
tkf1=zeros(Ny1,Nx1); % Old temperature fluid, K
tk2=zeros(Ny1,Nx1); % New temperature, K
vxp=zeros(Ny1,Nx1); % Solid Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Solid Vy in pressure nodes, m/s
vxpf=zeros(Ny1,Nx1); % Fluid Vx in pressure nodes, m/s
vypf=zeros(Ny1,Nx1); % Fluid Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Total Pressure, Pa
pf=zeros(Ny1,Nx1); % Fluid Pressure, Pa
PT0=zeros(Ny1,Nx1); % Old total Pressure, Pa
PF0=zeros(Ny1,Nx1); % Old fluid Pressure, Pa
ETAPHI=zeros(Ny1,Nx1); % Bulk Viscosity, Pa*s
PHI=zeros(Ny1,Nx1); % porosity
APHI=zeros(Ny1,Nx1); % Dln((1-PHI)/PHI)/Dt
XBF=zeros(Ny1,Nx1); % Basalt fraction in the fluid
XBS=zeros(Ny1,Nx1); % Basalt fraction in the solid
XHF=zeros(Ny1,Nx1); % Water fraction in the fluid
XHS=zeros(Ny1,Nx1); % Water fraction in the solid
XHT=zeros(Ny1,Nx1); % Total water fraction
DIVH=zeros(Ny1,Nx1); % Divergence of water flux

% Define markers
Nxmc=4; % Number of markers per cell in the horizontal direction
Nymc=4; % Number of markers per cell in the vertical direction
Nxm=(Nx-1)*Nxmc; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*Nymc; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
tkm=zeros(1,marknum); % Marker temperature, K
dtkfsm=zeros(1,marknum); % Marker fluid-solid temperature difference, K
rhocptotalm=zeros(1,marknum); % Marker volumetric Cp, J/m^3/K
phim=zeros(1,marknum); % Marker porosity
phicomm=zeros(1,marknum); % Marker compositional porosity
phinewm=zeros(1,marknum); % New marker porosity
XBsolidm0=zeros(1,marknum); % Old basalt content in solid
XHsolidm0=zeros(1,marknum); % Old water content in solid
XHsolidm=zeros(1,marknum); % New water content in solid
XHtotalm=zeros(1,marknum); % Total water content
pfm0=zeros(1,marknum); % Old fluid pressure, Pa
dxbf0=zeros(1,marknum); 
dxhf0=zeros(1,marknum); 
etasolidm=zeros(1,marknum); % Solid viscosity, Pa*s

% Define properties of materials:
% Shear Modulus, Pa
mu=8e+10; %80 GPa
% Burgers vector, m
b=5e-10; % 0.5 nm
% Grain size, m
h=1e-3; % 1 mm
% Dislocation creep (flowlaw) : [1=plagiocase 2=diabase 3=maficgranulite]
flow_law=1;
A_dis= [3.3e-4	2.0e-4  1.4e+4]; % Preexponential factor, MPa^-n*S^-n 	
n_dis= [3.2	    3.4     4.2   ]; % stress exponent
Ea_dis=[238000	260000  445000]; % Activation energy, J/mol
Va_dis=[0  	    0       0     ]; % Activation volume, J/bar/mol
% Diffusion creep		
A_dif= [8.7e+15	5.3e+15]; %1/s
m_dif= [-2.5	-2.5   ]; % grain size exponent
Ea_dif=[300000	240000 ]; % J/mol
Va_dif=[6e-6	5e-6   ]; % m^3/mol

%             mantle  plume  lith.   air
rhm         = [2      2      2       3000  ]; % Thermal conductivity, W/m/K
ksolidm     = [3      3      3       3000  ]; % Thermal conductivity, W/m/K
hrsolidm    = [2e-8   3e-8   2e-8    0     ]; % Radiogenic heat production, W/m^3
kphim0      = [1e-13  1e-13  1e-13   1e-13 ]; % Standard permeability, m^2

% Fluid properties (water)
%etafluid = 1e-5;     % Pa·s, viscosity of water
etafluid = 1e-3;     % Pa·s, viscosity of water
%etafluid = 1e1;     % Pa·s, viscosity of water
hrfluid  = 0;        % W/m^3, no heat production
kfluid   = 0.6;      % W/m/K, thermal conductivity of water

phim0=0.001; % standard porosity
phimin=1e-4; % Min porosity
phimax=1-phimin; % Max porosity
etamin=1e+16; % Min viscosity
etamax=1e+24; % Max viscosity
kphimin=1e-20; % Min permeability
kphimax=1e-10; % Max permeability

% Thermodynamic model: calibrated in this study
% Molar masses of oxides: g/mol
%       SiO2   MgO   FeO   Fe2O3   CaO   Al2O3   TiO2   MnO   K20   Na2O    H2O
Mi  =  [60.08 40.304 71.85  159.7  56.08 101.96  79.87  86.94 94.20 61.98  18.015]';

% Wt% Compositions of components: B=enriched ("basalt"), L=depleted("lherzolite"), H=water
% Lherzolite (Dummy)
X_L_wt=[ 0      0     0      0      0      0       0     0     0      0       0];
% Basalt (data from different missions; 1-6)
X_Venus_wt = [ ...
        50.8   12.8   7.7    0.0   8.3   17.8   0.22  0.16  0.11  2.0   0.0;  % 1: Kargel et al. (1993), VEGA
        45.6   11.5   7.7    0.0   7.5   16.0   0.20  0.14  0.10  2.0   0.0;  % 2: Surkov et al. (1986), VEGA-2
        48.7    8.1   7.92   0.98 10.3   17.9   1.30  8.10  0.20  2.4   0.0;  % 3: Semprich et al. (2024), Venusian Basalt
        45.1   11.4   9.3    0.0   7.1   15.8   1.60  0.20  4.0   2.0   0.0;  % 4: Filiberto et al. (2014); Fegley (2003), Venera-13
        48.7    8.1   8.8    0.0  10.3   17.9   1.30  0.20  0.2   2.4   0.0;  % 5: Filiberto et al. (2014); Fegley (2003), Venera-14
        45.0   37.8   8.05   0.0   3.55  4.45  0.201  0.135 0.029 0.36  0.0   % 6: McDonough & Sun (1995), BSE
];
X_B_wt=X_Venus_wt(Venus_composition, :);
% Water
X_H_wt=[ 0      0     0      0      0      0       0     0     0      0     100];


% Convert to mole fractions
% Basalt
X_B_mass = X_B_wt / 100;
n_B      = X_B_mass ./ Mi;
X_B_mol  = n_B / sum(n_B);
% Lherzolite
X_L_mass = X_L_wt / 100;
if sum(X_L_mass) > 0
    n_L     = X_L_mass ./ Mi;
    X_L_mol = n_L / sum(n_L);
else
    X_L_mol = zeros(size(Mi));    
end
% Water component
X_H_mass = X_H_wt / 100; 
n_H      = X_H_mass ./ Mi;
X_H_mol  = n_H / sum(n_H);       


idx_ox   = 1:10;                % oxides
Mi_ox    = Mi(idx_ox);          % molar masses of oxides
X_B_mol_ox = X_B_mol(idx_ox);   % mole fractions of oxides

% Renormalize to "dry" mole fractions
X_dry_mol = X_B_mol_ox / sum(X_B_mol_ox);

% Average molar mass of dry solid [kg/mol]
Mdry = sum(X_dry_mol .* Mi_ox);
% Convert molar masses to kg/mol
Mdry = Mdry / 1000;                % kg/mol
Mf   = 18.01528 / 1000;            % kg/mol  (water molar mass)


% Thermal boundary conditions
tktop=740;
tkbottom=840;
% Pressure boundary condition
ptop=95.6e+5; % Venus surface pressure, Pa
% Define marker coordinates, temperature and material type
an_width_x      = 30e3;         % full width in x: 30 km
an_width_y      = 20e3;         % full height in y: 10 km
an_half_x       = an_width_x/2; % 15 km
an_half_y       = an_width_y/2; % 5 km
% Gaussian standard deviations 
an_sigma_x      = an_half_x/2;  % ~7.5 km
an_sigma_y      = an_half_y/2;  % ~2.5 km
% Center of anomaly: horizontally centered, near the bottom
an_x_center     = xsize/2;              % center in x
an_y_center     = ysize;    % a bit above bottom, fully inside box
% Amplitudes of anomaly
an_dT_gauss     = 300.0;    % +300 K temperature anomaly
an_dXH_gauss    = 0.03;        % +0.05 in XHtotalm

m = 1; % Marker counter
for jm = 1:1:Nxm
    for im = 1:1:Nym
        % Define marker coordinates
        xm(m) = dxm/2 + (jm-1)*dxm + (rand-0.5)*dxm;
        ym(m) = dym/2 + (im-1)*dym + (rand-0.5)*dym;

        if (ym(m) < 5000)
            % Sticky air
            tm(m)        = 4;          % Material type
            tkm(m)       = tktop;        % Temperature
            phim(m)      = phimin;     % Porosity
            XHtotalm(m)  = 0;          % Water
            etasolidm(m) = 1e17;       % Initial viscosity in sticky air
        else
            tm(m)        = 1;
            tkm(m)       = tktop + ((ym(m)-5000)/20000)*100;
            phim(m)      = phim0;
            XHtotalm(m)  = 0.09;
            etasolidm(m) = 1e21;
        end
                
        
        if (ym(m) >= 5000)
            an_dx = xm(m) - an_x_center;
            an_dy = ym(m) - an_y_center;
            % Elliptical Gaussian
            an_r2 = (an_dx/an_sigma_x)^2 + (an_dy/an_sigma_y)^2;
            an_G  = exp(-0.5 * an_r2);
            if(an_G > 0.15)
                tm(m) = 2; % Plume material
                tkm(m)      = tkm(m)      + an_dT_gauss;
                % Water anomaly
                XHtotalm(m) = XHtotalm(m) + an_dXH_gauss;
            end
        end

        % Set solid hydration with interpolation from database and inital total hydration
        % approximate pressure for database
        if ym(m) <= 5000 
            ptnm  = ptop;
        else
            depth  = ym(m) - 5000;
            ptnm   = ptop + 3100 * gy * depth;
        end

        % Thermodynamic comptations
        T=tkm(m); % Temperature, K
        PB=(ptnm)*1e-9; % Pressure
        Wth20=XHtotalm(m); % Water content, mol %

        nt=fix((T-tkmin)/tkstp)+1;
        np=fix((PB-pbmin)/pbstp)+1;
        nwt=fix((Wth20-wth20min)/Wth20stp)+1;

        if(nt<1) 
            nt=1; 
        elseif(nt>tknum-1) 
            nt=tknum-1;
        end
        if(np<1) 
            np=1;
        elseif(np>pbnum-1) 
            np=pbnum-1;
        end
        if(nwt<1) 
            nwt=1;
        elseif(nwt>nwtnum-1) 
            nwt=nwtnum-1;
        end
        % Compute relative distances from the marker to the P-T-node
        dtt=(T-tkj(nt))/tkstp;
        dp=(PB-pbi(np))/pbstp;
        dwt=(Wth20-wth20k(nwt))/Wth20stp;

        c000 = XH2Osolid(np,  nt,   nwt  );
        c100 = XH2Osolid(np+1,nt,   nwt  );
        c010 = XH2Osolid(np,  nt+1, nwt  );
        c110 = XH2Osolid(np+1,nt+1, nwt  );
        c001 = XH2Osolid(np,  nt,   nwt+1);
        c101 = XH2Osolid(np+1,nt,   nwt+1);
        c011 = XH2Osolid(np,  nt+1, nwt+1);
        c111 = XH2Osolid(np+1,nt+1, nwt+1);
        % Corner weights
        w000 = (1-dp)*(1-dtt)*(1-dwt);
        w100 =    dp *(1-dtt)*(1-dwt);
        w010 = (1-dp)*   dtt *(1-dwt);
        w110 =    dp *   dtt *(1-dwt);
        w001 = (1-dp)*(1-dtt)*   dwt;
        w101 =    dp *(1-dtt)*   dwt;
        w011 = (1-dp)*   dtt *   dwt;
        w111 =    dp *   dtt *   dwt;
        % Hydration of solid (molar frac)
        XHsolidm0(m) = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                    w001*c001 + w101*c101 + w011*c011 + w111*c111;

        % Pyrolithe (no pyrolithe in this model)
        %XBsolidm0(m)=0.167964074*(1-XHsolidm0(m));  % pirolite
        XBsolidm0(m)=(1-XHsolidm0(m));
        
        % Update marker counter
        m = m + 1;
    end
end

XBsolidm=XBsolidm0;
XHsolidm=XHsolidm0;
phicomm=phim;
phinewm=phim;

% Introducing scaled pressure
pscale=1e+21/dx;

% Define global matrixes 
% Hydro-Mechanical solution: L(), R()
N=Nx1*Ny1*6; % Global number of unknowns
LL=zeros(N*44,1); % Line index in L
CL=zeros(N*44,1); % Column index in L 
VL=zeros(N*44,1); % Coefficient in L 
R=zeros(N,1); % Vector of the right parts of equations
L=sparse(N,N); % Matrix of coefficients (left part)
% Thermal solution: LT(), RT()
N=Nx1*Ny1; % Global number of unknowns
LLT=zeros(N*5,1); % Line index in L
CLT=zeros(N*5,1); % Column index in L 
VLT=zeros(N*5,1); % Coefficient in L 
RT=zeros(N,1); % Vector of right parts
LT=sparse(N,N); % Matrix of coefficients (left part)

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;
% Hydraulic boundary conditions: free slip=-1; No Slip=1
bcfleft=-1;
bcfright=-1;
bcftop=-1;
bcfbottom=-1;

% Timestepping
dt=1e+10; % initial timestep, s
timesum=0; % Time counter, s
dtkoef=1.2; % timestep increment
dxymax=0.1; % max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DTmax=20; % max temperature change per time step, K
dsubgridt=1; % subgrid diffusion parameter
dphimax=0.01; % max porosity ratio change per time step
pferrmax=10000; % Error limit to exit eterations, Pa
pfkoef=0.5; % Koeeficient for using pressure from the previous iteration
titermax=200; % Maximal number of global iterations
titermax1=20; % Maximal initial number of global iterations
titerstep1=5; % Number of timesteps for initial number of global iterations 
timestep=1;
savefigstep=50;  % .jpeg storage periodicity
savematstep=50;  % .mat storage periodicity
% pause(100)
end
desktopyn=1; % Figures y(1)/n(0)
for timestep=timestep:1:1000000

% Save old Ptotal for HA
PT00=PT0;

% HTM iteration
for titer=1:1:titermax

% Interpolate properties from markers to nodes
% Basic nodes
ETASUM=zeros(Ny,Nx);
TYPSUM=zeros(Ny,Nx);
%XBTSUM=zeros(Ny,Nx);
XHTSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
RHOFXSUM=zeros(Ny1,Nx1);
KXSUM=zeros(Ny1,Nx1);
PHIXSUM=zeros(Ny1,Nx1);
RXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
RHOFYSUM=zeros(Ny1,Nx1);
KYSUM=zeros(Ny1,Nx1);
PHIYSUM=zeros(Ny1,Nx1);
RYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
% P-Nodes
ETAPSUM=zeros(Ny1,Nx1);
ETASPSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
ALPHASUM=zeros(Ny1,Nx1);
ALPHASSUM=zeros(Ny1,Nx1);
BETTASSUM=zeros(Ny1,Nx1);
ALPHAFSUM=zeros(Ny1,Nx1);
BETTAFSUM=zeros(Ny1,Nx1);
HRSUM=zeros(Ny1,Nx1);
TKSUM=zeros(Ny1,Nx1);
DTKSUM=zeros(Ny1,Nx1);
PHISUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
XHSSUM=zeros(Ny1,Nx1);
XHTSUM=zeros(Ny1,Nx1);
RHOWATERSUM=zeros(Ny1,Nx1);
RHOSOLIDSUM=zeros(Ny1,Nx1);


for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
    
    % Compute marker parameters
    % Compute density of solid and fluid
    % Interpolate fluid pressure and temperature to the marker
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;

    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Interpolate nodal total pressure and temperature to the marker
    ptnm=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
    tknm=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    ptnm=max(ptnm,ptop);
    if(timestep==1 && titer==1)
        % approximate pressure for database (only for the first timestep)
        if (ym(m) <= 5000) 
            ptnm  = ptop;
        else
            depth  = ym(m) - 5000;
            ptnm   = ptop + 3100 * gy * depth;
        end
    end
    xhtnm=XHtotalm(m);
    Wth20=xhtnm;

    % Thermodynamic computations
    T=tknm; % Temperature, K
    if(T<=0)
        T=tkm(m);
    end
    PB=(ptnm)*1e-9; % Pressure, Pa=1e+9 GPa

    nt=fix((T-tkmin)/tkstp)+1;
    np=fix((PB-pbmin)/pbstp)+1;
    nwt=fix((Wth20-wth20min)/Wth20stp)+1;

    if(nt<1) 
        nt=1; 
    elseif(nt>tknum-1) 
        nt=tknum-1;
    end
    if(np<1) 
        np=1;
    elseif(np>pbnum-1) 
        np=pbnum-1;
    end
    if(nwt<1) 
        nwt=1;
    elseif(nwt>nwtnum-1) 
        nwt=nwtnum-1;
    end
    % Compute relative distances from the marker to the P-T-node
    dtt=(T-tkj(nt))/tkstp;
    dp=(PB-pbi(np))/pbstp;
    dwt=(Wth20-wth20k(nwt))/Wth20stp;

    % Compute weights
    w000 = (1-dp)*(1-dtt)*(1-dwt);
    w100 =    dp *(1-dtt)*(1-dwt);
    w010 = (1-dp)*   dtt *(1-dwt);
    w110 =    dp *   dtt *(1-dwt);
    w001 = (1-dp)*(1-dtt)*   dwt;
    w101 =    dp *(1-dtt)*   dwt;
    w011 = (1-dp)*   dtt *   dwt;
    w111 =    dp *   dtt *   dwt;


    % Water density, alpha, betta interpolation: Gerya et al., 2004, PCM
    % Water density, alpha, betta interpolation: Z. Jing, S. Karato / Geochimica et Cosmochimica Acta 85 (2012) 357–372 
    % Corner densities
    c000 = rholiquid(np,  nt,   nwt  );
    c100 = rholiquid(np+1,nt,   nwt  );
    c010 = rholiquid(np,  nt+1, nwt  );
    c110 = rholiquid(np+1,nt+1, nwt  );
    c001 = rholiquid(np,  nt,   nwt+1);
    c101 = rholiquid(np+1,nt,   nwt+1);
    c011 = rholiquid(np,  nt+1, nwt+1);
    c111 = rholiquid(np+1,nt+1, nwt+1);
    % Density of water
    rho_water = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
           w001*c001 + w101*c101 + w011*c011 + w111*c111;
    % Corner alpha
    c000 = a_liquid(np,  nt,   nwt  );
    c100 = a_liquid(np+1,nt,   nwt  );
    c010 = a_liquid(np,  nt+1, nwt  );
    c110 = a_liquid(np+1,nt+1, nwt  );
    c001 = a_liquid(np,  nt,   nwt+1);
    c101 = a_liquid(np+1,nt,   nwt+1);
    c011 = a_liquid(np,  nt+1, nwt+1);
    c111 = a_liquid(np+1,nt+1, nwt+1);
    % Thermal expansion
    alpha_water = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                  w001*c001 + w101*c101 + w011*c011 + w111*c111;
    % Corner betta
    c000 = b_liquid(np,  nt,   nwt  );
    c100 = b_liquid(np+1,nt,   nwt  );
    c010 = b_liquid(np,  nt+1, nwt  );
    c110 = b_liquid(np+1,nt+1, nwt  );
    c001 = b_liquid(np,  nt,   nwt+1);
    c101 = b_liquid(np+1,nt,   nwt+1);
    c011 = b_liquid(np,  nt+1, nwt+1);
    c111 = b_liquid(np+1,nt+1, nwt+1);
    % Thermal compressibility
    betta_water = (w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                   w001*c001 + w101*c101 + w011*c011 + w111*c111)*1e-9;
    
    % Compute wet solid density, alpha, betta
    % Corner densities
    c000 = rhosolid(np,  nt,   nwt  );
    c100 = rhosolid(np+1,nt,   nwt  );
    c010 = rhosolid(np,  nt+1, nwt  );
    c110 = rhosolid(np+1,nt+1, nwt  );
    c001 = rhosolid(np,  nt,   nwt+1);
    c101 = rhosolid(np+1,nt,   nwt+1);
    c011 = rhosolid(np,  nt+1, nwt+1);
    c111 = rhosolid(np+1,nt+1, nwt+1);
    % Density of solid
    rho_solid = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                   w001*c001 + w101*c101 + w011*c011 + w111*c111;
    % Corner alpha
    c000 = a_solid(np,  nt,   nwt  );
    c100 = a_solid(np+1,nt,   nwt  );
    c010 = a_solid(np,  nt+1, nwt  );
    c110 = a_solid(np+1,nt+1, nwt  );
    c001 = a_solid(np,  nt,   nwt+1);
    c101 = a_solid(np+1,nt,   nwt+1);
    c011 = a_solid(np,  nt+1, nwt+1);
    c111 = a_solid(np+1,nt+1, nwt+1);
    % Thermal expansion
    alpha_solid = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                     w001*c001 + w101*c101 + w011*c011 + w111*c111;
    % Corner betta
    c000 = b_solid(np,  nt,   nwt  );
    c100 = b_solid(np+1,nt,   nwt  );
    c010 = b_solid(np,  nt+1, nwt  );
    c110 = b_solid(np+1,nt+1, nwt  );
    c001 = b_solid(np,  nt,   nwt+1);
    c101 = b_solid(np+1,nt,   nwt+1);
    c011 = b_solid(np,  nt+1, nwt+1);
    c111 = b_solid(np+1,nt+1, nwt+1);
    % Thermal compressibilty
    betta_solid = (w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                     w001*c001 + w101*c101 + w011*c011 + w111*c111) * 1e-9;

    % System density
    %From database
    c000 = rhosys(np,  nt,   nwt  );
    c100 = rhosys(np+1,nt,   nwt  );
    c010 = rhosys(np,  nt+1, nwt  );
    c110 = rhosys(np+1,nt+1, nwt  );
    c001 = rhosys(np,  nt,   nwt+1);
    c101 = rhosys(np+1,nt,   nwt+1);
    c011 = rhosys(np,  nt+1, nwt+1);
    c111 = rhosys(np+1,nt+1, nwt+1);
    % Density of solid
    rho_sys = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
              w001*c001 + w101*c101 + w011*c011 + w111*c111;

    % Cptotal
    c000 = cptot(np,  nt,   nwt  );
    c100 = cptot(np+1,nt,   nwt  );
    c010 = cptot(np,  nt+1, nwt  );
    c110 = cptot(np+1,nt+1, nwt  );
    c001 = cptot(np,  nt,   nwt+1);
    c101 = cptot(np+1,nt,   nwt+1);
    c011 = cptot(np,  nt+1, nwt+1);
    c111 = cptot(np+1,nt+1, nwt+1);
    % Total heat capacity of solid
    cptot_sys = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
              w001*c001 + w101*c101 + w011*c011 + w111*c111;
    cptot_sys = cptot_sys * 1000; % J/kg/K

    % Porosity at equilibrium
    % Density of solid and fluid
    rhosolidm0 = rho_solid;  %Hydrated solid density
    rhofluidm0 = rho_water;  %Water density
    rhototalm = rho_sys;     %System density
    
    denom_phi = (rhofluidm0 - rhosolidm0);
    if abs(denom_phi) > 0
        phi_eq = (rhototalm - rhosolidm0) / denom_phi;
        phi_eq = max(0, min(phi_eq, 1));
    else
        phi_eq = phicomm(m);
    end

    phicomm(m) = phi_eq;
    phim(m) = max(phimin, min(phicomm(m), phimax));

    % Compute properties
    kphim    = min(max(kphim0(tm(m))*(phim(m)/phim0)^3/((1-phim(m))/(1-phim0))^2,kphimin),kphimax);
    % bulk Cp and alpha (solid+water)
    rhocptotalm(m) = cptot_sys * rho_sys;
    alphatotalm    = (1-phim(m))*alpha_solid + phim(m)*alpha_water;
    % Viscosity
    etatotalm=min(max(etasolidm(m)*exp(-28*phim(m)),etamin),etamax);

    hrtotalm=hrsolidm(tm(m))*(1-phim(m))+hrfluid*phim(m);
    ktotalm=(ksolidm(tm(m))*kfluid/2+((ksolidm(tm(m))*(3*phim(m)-2)+...
        kfluid*(1-3*phim(m)))^2)/16)^0.5-(ksolidm(tm(m))*(3*phim(m)-2)+...
        kfluid*(1-3*phim(m)))/4;
    % pause(100)

if(tm(m)==4)
    % Sticky Air
    hrtotalm=hrsolidm(tm(m));
    ktotalm=ksolidm(tm(m));
    rhototalm=100;
    rhocptotalm(m)=3e+6;
    rho_water=0;
    rho_solid=0;
    etatotalm=etasolidm(m);
    kphim=1e-17; %Permeability
    XHtotalm(m)=XHsolidm0(m);
    alphatotalm=0;
    alpha_solid=0;
    alpha_water=0;
    betta_solid=1e-10;
    betta_water=1e-10;
end


    % Interpolation to basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    ETASUM(i,j)=ETASUM(i,j)+etatotalm*wtmij;
    TYPSUM(i,j)=TYPSUM(i,j)+tm(m)*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETASUM(i+1,j)=ETASUM(i+1,j)+etatotalm*wtmi1j;
    TYPSUM(i+1,j)=TYPSUM(i+1,j)+tm(m)*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETASUM(i,j+1)=ETASUM(i,j+1)+etatotalm*wtmij1;
    TYPSUM(i,j+1)=TYPSUM(i,j+1)+tm(m)*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etatotalm*wtmi1j1;
    TYPSUM(i+1,j+1)=TYPSUM(i+1,j+1)+tm(m)*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;    
    
    
    % Interpolation to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhototalm*wtmij;
    RHOFXSUM(i,j)=RHOFXSUM(i,j)+rhofluidm0*wtmij;
    KXSUM(i,j)=KXSUM(i,j)+ktotalm*wtmij;
    PHIXSUM(i,j)=PHIXSUM(i,j)+phim(m)*wtmij;
    RXSUM(i,j)=RXSUM(i,j)+etafluid/kphim*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFXSUM(i+1,j)=RHOFXSUM(i+1,j)+rhofluidm0*wtmi1j;
    KXSUM(i+1,j)=KXSUM(i+1,j)+ktotalm*wtmi1j;
    PHIXSUM(i+1,j)=PHIXSUM(i+1,j)+phim(m)*wtmi1j;
    RXSUM(i+1,j)=RXSUM(i+1,j)+etafluid/kphim*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhototalm*wtmij1;
    RHOFXSUM(i,j+1)=RHOFXSUM(i,j+1)+rhofluidm0*wtmij1;
    KXSUM(i,j+1)=KXSUM(i,j+1)+ktotalm*wtmij1;
    PHIXSUM(i,j+1)=PHIXSUM(i,j+1)+phim(m)*wtmij1;
    RXSUM(i,j+1)=RXSUM(i,j+1)+etafluid/kphim*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFXSUM(i+1,j+1)=RHOFXSUM(i+1,j+1)+rhofluidm0*wtmi1j1;
    KXSUM(i+1,j+1)=KXSUM(i+1,j+1)+ktotalm*wtmi1j1;
    RXSUM(i+1,j+1)=RXSUM(i+1,j+1)+etafluid/kphim*wtmi1j1;
    PHIXSUM(i+1,j+1)=PHIXSUM(i+1,j+1)+phim(m)*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhototalm*wtmij;
    RHOFYSUM(i,j)=RHOFYSUM(i,j)+rhofluidm0*wtmij;
    KYSUM(i,j)=KYSUM(i,j)+ktotalm*wtmij;
    PHIYSUM(i,j)=PHIYSUM(i,j)+phim(m)*wtmij;
    RYSUM(i,j)=RYSUM(i,j)+etafluid/kphim*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFYSUM(i+1,j)=RHOFYSUM(i+1,j)+rhofluidm0*wtmi1j;
    KYSUM(i+1,j)=KYSUM(i+1,j)+ktotalm*wtmi1j;
    PHIYSUM(i+1,j)=PHIYSUM(i+1,j)+phim(m)*wtmi1j;
    RYSUM(i+1,j)=RYSUM(i+1,j)+etafluid/kphim*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhototalm*wtmij1;
    RHOFYSUM(i,j+1)=RHOFYSUM(i,j+1)+rhofluidm0*wtmij1;
    KYSUM(i,j+1)=KYSUM(i,j+1)+ktotalm*wtmij1;
    PHIYSUM(i,j+1)=PHIYSUM(i,j+1)+phim(m)*wtmij1;
    RYSUM(i,j+1)=RYSUM(i,j+1)+etafluid/kphim*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFYSUM(i+1,j+1)=RHOFYSUM(i+1,j+1)+rhofluidm0*wtmi1j1;
    KYSUM(i+1,j+1)=KYSUM(i+1,j+1)+ktotalm*wtmi1j1;
    PHIYSUM(i+1,j+1)=PHIYSUM(i+1,j+1)+phim(m)*wtmi1j1;
    RYSUM(i+1,j+1)=RYSUM(i+1,j+1)+etafluid/kphim*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to P-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    XHSSUM(i,j)=XHSSUM(i,j)+XHsolidm0(m)*wtmij;
    XHTSUM(i,j)=XHTSUM(i,j)+XHtotalm(m)*wtmij;
    ETAPSUM(i,j)=ETAPSUM(i,j)+etatotalm*wtmij;
    ETASPSUM(i,j)=ETASPSUM(i,j)+etasolidm(m)*wtmij;
    RHOSUM(i,j)=RHOSUM(i,j)+rhototalm*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocptotalm(m)*wtmij;
    ALPHASUM(i,j)=ALPHASUM(i,j)+alphatotalm*wtmij;
    ALPHASSUM(i,j)=ALPHASSUM(i,j)+alpha_solid*wtmij;
    ALPHAFSUM(i,j)=ALPHAFSUM(i,j)+alpha_water*wtmij;
    BETTASSUM(i,j)=BETTASSUM(i,j)+betta_solid*wtmij;
    BETTAFSUM(i,j)=BETTAFSUM(i,j)+betta_water*wtmij;
    HRSUM(i,j)=HRSUM(i,j)+hrtotalm*wtmij;
    TKSUM(i,j)=TKSUM(i,j)+tkm(m)*rhocptotalm(m)*wtmij;
    DTKSUM(i,j)=DTKSUM(i,j)+dtkfsm(m)*rhocptotalm(m)*wtmij;
    PHISUM(i,j)=PHISUM(i,j)+phim(m)*wtmij;
    RHOWATERSUM(i,j)=RHOWATERSUM(i,j)+rho_water*wtmij;
    RHOSOLIDSUM(i,j)=RHOSOLIDSUM(i,j)+rho_solid*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    XHSSUM(i+1,j)=XHSSUM(i+1,j)+XHsolidm0(m)*wtmi1j;
    XHTSUM(i+1,j)=XHTSUM(i+1,j)+XHtotalm(m)*wtmi1j;
    ETAPSUM(i+1,j)=ETAPSUM(i+1,j)+etatotalm*wtmi1j;
    ETASPSUM(i+1,j)=ETASPSUM(i+1,j)+etasolidm(m)*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhototalm*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocptotalm(m)*wtmi1j;
    ALPHASUM(i+1,j)=ALPHASUM(i+1,j)+alphatotalm*wtmi1j;
    ALPHASSUM(i+1,j)=ALPHASSUM(i+1,j)+alpha_solid*wtmi1j;
    ALPHAFSUM(i+1,j)=ALPHAFSUM(i+1,j)+alpha_water*wtmi1j;
    BETTASSUM(i+1,j)=BETTASSUM(i+1,j)+betta_solid*wtmi1j;
    BETTAFSUM(i+1,j)=BETTAFSUM(i+1,j)+betta_water*wtmi1j;
    HRSUM(i+1,j)=HRSUM(i+1,j)+hrtotalm*wtmi1j;
    TKSUM(i+1,j)=TKSUM(i+1,j)+tkm(m)*rhocptotalm(m)*wtmi1j;
    DTKSUM(i+1,j)=DTKSUM(i+1,j)+dtkfsm(m)*rhocptotalm(m)*wtmi1j;
    PHISUM(i+1,j)=PHISUM(i+1,j)+phim(m)*wtmi1j;
    RHOWATERSUM(i+1,j)=RHOWATERSUM(i+1,j)+rho_water*wtmi1j;
    RHOSOLIDSUM(i+1,j)=RHOSOLIDSUM(i+1,j)+rho_solid*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    XHSSUM(i,j+1)=XHSSUM(i,j+1)+XHsolidm0(m)*wtmij1;
    XHTSUM(i,j+1)=XHTSUM(i,j+1)+XHtotalm(m)*wtmij1;
    ETAPSUM(i,j+1)=ETAPSUM(i,j+1)+etatotalm*wtmij1;
    ETASPSUM(i,j+1)=ETASPSUM(i,j+1)+etasolidm(m)*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhototalm*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocptotalm(m)*wtmij1;
    ALPHASUM(i,j+1)=ALPHASUM(i,j+1)+alphatotalm*wtmij1;
    ALPHASSUM(i,j+1)=ALPHASSUM(i,j+1)+alpha_solid*wtmij1;
    ALPHAFSUM(i,j+1)=ALPHAFSUM(i,j+1)+alpha_water*wtmij1;
    BETTASSUM(i,j+1)=BETTASSUM(i,j+1)+betta_solid*wtmij1;
    BETTAFSUM(i,j+1)=BETTAFSUM(i,j+1)+betta_water*wtmij1;
    HRSUM(i,j+1)=HRSUM(i,j+1)+hrtotalm*wtmij1;
    TKSUM(i,j+1)=TKSUM(i,j+1)+tkm(m)*rhocptotalm(m)*wtmij1;
    DTKSUM(i,j+1)=DTKSUM(i,j+1)+dtkfsm(m)*rhocptotalm(m)*wtmij1;
    PHISUM(i,j+1)=PHISUM(i,j+1)+phim(m)*wtmij1;
    RHOWATERSUM(i,j+1)=RHOWATERSUM(i,j+1)+rho_water*wtmij1;
    RHOSOLIDSUM(i,j+1)=RHOSOLIDSUM(i,j+1)+rho_solid*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    XHSSUM(i+1,j+1)=XHSSUM(i+1,j+1)+XHsolidm0(m)*wtmi1j1;
    XHTSUM(i+1,j+1)=XHTSUM(i+1,j+1)+XHtotalm(m)*wtmi1j1;
    ETAPSUM(i+1,j+1)=ETAPSUM(i+1,j+1)+etatotalm*wtmi1j1;
    ETASPSUM(i+1,j+1)=ETASPSUM(i+1,j+1)+etasolidm(m)*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocptotalm(m)*wtmi1j1;
    ALPHASUM(i+1,j+1)=ALPHASUM(i+1,j+1)+alphatotalm*wtmi1j1;
    ALPHASSUM(i+1,j+1)=ALPHASSUM(i+1,j+1)+alpha_solid*wtmi1j1;
    ALPHAFSUM(i+1,j+1)=ALPHAFSUM(i+1,j+1)+alpha_water*wtmi1j1;
    BETTASSUM(i+1,j+1)=BETTASSUM(i+1,j+1)+betta_solid*wtmi1j1;
    BETTAFSUM(i+1,j+1)=BETTAFSUM(i+1,j+1)+betta_water*wtmi1j1;
    HRSUM(i+1,j+1)=HRSUM(i+1,j+1)+hrtotalm*wtmi1j1;
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+tkm(m)*rhocptotalm(m)*wtmi1j1;
    DTKSUM(i+1,j+1)=DTKSUM(i+1,j+1)+dtkfsm(m)*rhocptotalm(m)*wtmi1j1;
    PHISUM(i+1,j+1)=PHISUM(i+1,j+1)+phim(m)*wtmi1j1;
    RHOWATERSUM(i+1,j+1)=RHOWATERSUM(i+1,j+1)+rho_water*wtmi1j1;
    RHOSOLIDSUM(i+1,j+1)=RHOSOLIDSUM(i+1,j+1)+rho_solid*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
end
% Compute physical properties
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
            TYP(i,j)=TYPSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
            RHOFX(i,j)=RHOFXSUM(i,j)/WTXSUM(i,j);
            KX(i,j)=KXSUM(i,j)/WTXSUM(i,j);
            PHIX(i,j)=PHIXSUM(i,j)/WTXSUM(i,j);
            RX(i,j)=RXSUM(i,j)/WTXSUM(i,j);
        end
%         RX(i,j)=min(max(RX(i,j),etafluid/kphimax),etafluid/kphimin);
    end
end
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
            RHOFY(i,j)=RHOFYSUM(i,j)/WTYSUM(i,j);
            KY(i,j)=KYSUM(i,j)/WTYSUM(i,j);
            PHIY(i,j)=PHIYSUM(i,j)/WTYSUM(i,j);
            RY(i,j)=RYSUM(i,j)/WTYSUM(i,j);
        end
%         RY(i,j)=min(max(RY(i,j),etafluid/kphimax),etafluid/kphimin);
    end
end
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            XHS(i,j)=XHSSUM(i,j)/WTPSUM(i,j);
            XHT(i,j)=XHTSUM(i,j)/WTPSUM(i,j);
            PHI(i,j)=PHISUM(i,j)/WTPSUM(i,j);
            ETAP(i,j)=ETAPSUM(i,j)/WTPSUM(i,j);
            ETAPHI(i,j)=PHI(i,j)/(ETASPSUM(i,j)/WTPSUM(i,j));
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            RHOCP(i,j)=RHOCPSUM(i,j)/WTPSUM(i,j);
            ALPHA(i,j)=ALPHASUM(i,j)/WTPSUM(i,j);
            ALPHASOLID(i,j)=ALPHASSUM(i,j)/WTPSUM(i,j);
            ALPHAFLUID(i,j)=ALPHAFSUM(i,j)/WTPSUM(i,j);
            BETTASOLID(i,j)=BETTASSUM(i,j)/WTPSUM(i,j);
            BETTAFLUID(i,j)=BETTAFSUM(i,j)/WTPSUM(i,j);
            HR(i,j)=HRSUM(i,j)/WTPSUM(i,j);
            tk1(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
            tkf1(i,j)=(TKSUM(i,j)+DTKSUM(i,j))/RHOCPSUM(i,j);
            RHOWATER(i,j)=RHOWATERSUM(i,j)/WTPSUM(i,j);
            RHOSOLID(i,j)=RHOSOLIDSUM(i,j)/WTPSUM(i,j);
        end
    end
end

%%
% Save initial temperature array
if(timestep==1 && titer==1)
    tk00=tk1;
end
% Make adiabatic temperature profile
if(timestep==0)
    for j=2:1:Nx
        dtk1=RHOY(1,j)*dy/2*gy*ALPHA(2,j)/RHOCP(2,j);
        tk1(2,j)=tk00(2,j)+dtk1;
        for i=3:1:Ny
            if(titer==1)
                dtk1=dtk1+RHOY(i-1,2)*dy*gy*(ALPHA(i-1,2)/RHOCP(i-1,2)+ALPHA(i,2)/RHOCP(i,2))/2;
            else
                dtk1=dtk1+(pr(i,2)-pr(i-1,2))*(ALPHA(i-1,2)/RHOCP(i-1,2)+ALPHA(i,2)/RHOCP(i,2))/2;
            end
            tk1(i,j)=tk00(i,j)+dtk1;
        end
        % Bottom temperature
        dtk1=RHOY(Ny,j)*dy/2*gy*ALPHA(Ny,j)/RHOCP(Ny,j);
        tk1(Ny1,j)=tk1(Ny,j)+dtk1;
    end
% Recalculate bottom temperature
tkbottom=mean(tk1(Ny1,2:Nx));
end

% Applying thermal boundary conditions for interpolated temperature
% Upper boundary 
tk1(1,2:Nx)=2*tktop-tk1(2,2:Nx); % Constant temperature
% Lower boundary 
tk1(Ny1,2:Nx)=2*tkbottom-tk1(Ny,2:Nx); % Constant temperature
% Left boundary
tk1(:,1)=tk1(:,2); % Insulating boundary
% Right boundary
tk1(:,Nx1)=tk1(:,Nx); % Insulating boundary
% Reset new temperature array
if(titer==1)
    tk2=tk1;
end

% Compute mass transfer rate
DMPSUM=zeros(Ny1,Nx1);
DHPSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
% Mass transfer rate here (not in this model)
% P-nodes
DMP=zeros(Ny1,Nx1);
DHP=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            DMP(i,j)=DMPSUM(i,j)/WTPSUM(i,j);
            DHP(i,j)=DHPSUM(i,j)/WTPSUM(i,j);
        end
    end
end
DMPmin=min(min(DMP))    
DMPmax=max(max(DMP))    
DHPmin=min(min(DHP))    
DHPmax=max(max(DHP))    
 

% Hydro-Mechanical Solution
pf0=pf;
% Elastic part y/n
gggbkoef=1;
if(timestep==1)
    gggbkoef=0;
end
% Reset global counter
GL=1;
% Composing global matrixes L(), R() for Stokes and continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % P fluid
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition
            % Left BC, Right BC, Ghost nodes: 1*vx=0
            if(j==1 || j==Nx || j==Nx1)
                LL(GL)=kvx; CL(GL)=kvx; VL(GL)=1; GL=GL+1;
                R(kvx)=0;
            % Top boundary: free slip/no slip
            elseif(i==1)
                LL(GL)=kvx; CL(GL)=kvx; VL(GL)=1; GL=GL+1;
                LL(GL)=kvx; CL(GL)=kvx+6; VL(GL)=bctop; GL=GL+1;
                R(kvx)=0;
            else
            % Bottom boundary: free slip/no slip
                LL(GL)=kvx; CL(GL)=kvx-6; VL(GL)=bcbottom; GL=GL+1;
                LL(GL)=kvx; CL(GL)=kvx; VL(GL)=1; GL=GL+1;
                R(kvx)=0;
            end
        else
        % Internal points: x-Stokes eq.
        % ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
        %            Vx2
        %             |
        %        Vy1  |  Vy3
        %             |
        %     Vx1-P1-Vx3-P2-Vx5
        %             |
        %        Vy2  |  Vy4
        %             |
        %            Vx4
        %
        % Viscosity points
        ETA1=ETA(i-1,j);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i,j+1);
        % Left part
        LL(GL)=kvx; CL(GL)=kvx-Ny1*6; VL(GL)=ETAP1/dx^2; GL=GL+1; % Vx1
        LL(GL)=kvx; CL(GL)=kvx-6; VL(GL)=ETA1/dy^2; GL=GL+1; % Vx2
        LL(GL)=kvx; CL(GL)=kvx; VL(GL)=-(ETAP1+ETAP2)/dx^2-(ETA1+ETA2)/dy^2; GL=GL+1; % Vx3
        LL(GL)=kvx; CL(GL)=kvx+6; VL(GL)=ETA2/dy^2; GL=GL+1; % Vx4
        LL(GL)=kvx; CL(GL)=kvx+Ny1*6; VL(GL)=ETAP2/dx^2; GL=GL+1; % Vx5
        LL(GL)=kvx; CL(GL)=kvy; VL(GL)=ETAP1/dx/dy-ETA2/dx/dy; GL=GL+1; % Vy2
        LL(GL)=kvx; CL(GL)=kvy+Ny1*6; VL(GL)=-ETAP2/dx/dy+ETA2/dx/dy; GL=GL+1; % Vy4
        LL(GL)=kvx; CL(GL)=kvy-6; VL(GL)=-ETAP1/dx/dy+ETA1/dx/dy; GL=GL+1; % Vy1
        LL(GL)=kvx; CL(GL)=kvy+Ny1*6-6; VL(GL)=ETAP2/dx/dy-ETA1/dx/dy; GL=GL+1; % Vy3
        LL(GL)=kvx; CL(GL)=kpm; VL(GL)=pscale/dx; GL=GL+1; % P1
        LL(GL)=kvx; CL(GL)=kpm+Ny1*6; VL(GL)=-pscale/dx; GL=GL+1; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*gx;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Top BC, Bottom BC, Ghost nodes: 1*vy=0
            if(i==1 && j>1 && j<Nx1)
                % %            Vy1
                % %             |
                % %      Vx1--Pt,Pf--Vx2
                % %             |
                % %            Vy2
                % %
                % % Left part
                % LL(GL)=kvy; CL(GL)=kvx+6-Ny1*6; VL(GL)=-1/dx; GL=GL+1; % Vx1
                % LL(GL)=kvy; CL(GL)=kvx+6; VL(GL)=1/dx; GL=GL+1; % Vx2
                % LL(GL)=kvy; CL(GL)=kvy; VL(GL)=-1/dy; GL=GL+1; % Vy1
                % LL(GL)=kvy; CL(GL)=kvy+6; VL(GL)=1/dy; GL=GL+1; % Vy2
                % % Right part
                % R(kvy)=0;
                LL(GL)=kvy; CL(GL)=kvy; VL(GL)=1; GL=GL+1;
                R(kvy)=0;
            elseif(i==Ny || i==Ny1)
                LL(GL)=kvy; CL(GL)=kvy; VL(GL)=1; GL=GL+1;
                R(kvy)=0;
            % Left boundary: free slip/no slip
            elseif(j==1)
                LL(GL)=kvy; CL(GL)=kvy; VL(GL)=1; GL=GL+1;
                LL(GL)=kvy; CL(GL)=kvy+6*Ny1; VL(GL)=bcleft; GL=GL+1;
                R(kvy)=0;
            else
            % Right boundary: free slip/no slip
                LL(GL)=kvy; CL(GL)=kvy-6*Ny1; VL(GL)=bcright; GL=GL+1;
                LL(GL)=kvy; CL(GL)=kvy; VL(GL)=1; GL=GL+1;
                R(kvy)=0;
            end
      
        else
        % Internal points: y-Stokes eq.
        % ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
        %            Vy2
        %             |
        %         Vx1 P1 Vx3
        %             |
        %     Vy1----Vy3----Vy5
        %             |
        %         Vx2 P2 Vx4
        %             |
        %            Vy4
        %
        % Viscosity points
        % Viscosity points
        ETA1=ETA(i,j-1);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i+1,j);
        % Left part
        LL(GL)=kvy; CL(GL)=kvy-Ny1*6; VL(GL)=ETA1/dx^2; GL=GL+1; % Vy1
        LL(GL)=kvy; CL(GL)=kvy-6; VL(GL)=ETAP1/dy^2; GL=GL+1; % Vy2
        LL(GL)=kvy; CL(GL)=kvy; VL(GL)=-(ETAP1+ETAP2)/dy^2-(ETA1+ETA2)/dx^2; GL=GL+1; % Vy3
        LL(GL)=kvy; CL(GL)=kvy+6; VL(GL)=ETAP2/dy^2; GL=GL+1; % Vy4
        LL(GL)=kvy; CL(GL)=kvy+Ny1*6; VL(GL)=ETA2/dx^2; GL=GL+1; % Vy5
        LL(GL)=kvy; CL(GL)=kvx; VL(GL)=ETAP1/dx/dy-ETA2/dx/dy; GL=GL+1; % Vx3
        LL(GL)=kvy; CL(GL)=kvx+6; VL(GL)=-ETAP2/dx/dy+ETA2/dx/dy; GL=GL+1; % Vx4
        LL(GL)=kvy; CL(GL)=kvx-Ny1*6; VL(GL)=-ETAP1/dx/dy+ETA1/dx/dy; GL=GL+1; % Vx1
        LL(GL)=kvy; CL(GL)=kvx+6-Ny1*6; VL(GL)=ETAP2/dx/dy-ETA1/dx/dy; GL=GL+1; % Vx2
        LL(GL)=kvy; CL(GL)=kpm; VL(GL)=pscale/dy; GL=GL+1; % P1
        LL(GL)=kvy; CL(GL)=kpm+6; VL(GL)=-pscale/dy; GL=GL+1; % P2
        % Right part
        R(kvy)=-RHOY(i,j)*gy;
        end
        
        % Ptotal equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            % Ghost nodes: 1*P=0
            LL(GL)=kpm; CL(GL)=kpm; VL(GL)=1*pscale; GL=GL+1;
            R(kpm)=0; % Right part
        % Real BC
        elseif(i==2 && j==2 && gggbkoef==0)
            LL(GL)=kpm; CL(GL)=kpm; VL(GL)=1*pscale; GL=GL+1;
            R(kpm)=ptop; % Right part
        else
        % Grid steps
        % Internal points: continuity eq.
        % dVx/dx+dVy/dy+(Ptotal-Pfluid)/ETHAphi/(1-PHI)+BETTAd*(DsPt/dt-Kbw*DfPf/dt)=DMP+ALPHAs*DsT/Dt
        %            Vy1
        %             |
        %      Vx1--Pt,Pf--Vx2
        %             |
        %            Vy2
        %
        % Pores compressibility
        BETTAPHI=0;
        % Drained compressibility
        BETTADRAINED=(BETTAPHI+BETTASOLID(i,j))/(1-PHI(i,j)); 
        % Biott-Willis koefficient
        KBW=1-BETTASOLID(i,j)/BETTADRAINED; 
        % Left part
        LL(GL)=kpm; CL(GL)=kvx-Ny1*6; VL(GL)=-1/dx; GL=GL+1; % Vx1
        LL(GL)=kpm; CL(GL)=kvx; VL(GL)=1/dx; GL=GL+1; % Vx2
        LL(GL)=kpm; CL(GL)=kvy-6; VL(GL)=-1/dy; GL=GL+1; % Vy1
        LL(GL)=kpm; CL(GL)=kvy; VL(GL)=1/dy; GL=GL+1; % Vy2
        LL(GL)=kpm; CL(GL)=kpm; VL(GL)=pscale*(ETAPHI(i,j)/(1-PHI(i,j))+gggbkoef*BETTADRAINED/dt); GL=GL+1; % Ptotal
        LL(GL)=kpm; CL(GL)=kpf; VL(GL)=-pscale*(ETAPHI(i,j)/(1-PHI(i,j))+gggbkoef*BETTADRAINED*KBW/dt); GL=GL+1; % Pfluid
        % Right part
        R(kpm)=DMP(i,j)+gggbkoef*(BETTADRAINED*(PT0(i,j)-KBW*PF0(i,j))/dt+ALPHASOLID(i,j)*(tk2(i,j)-tk1(i,j))/dt);
        end

        % qxDarcy equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Left BC, Right BC, Ghost nodes: 1*qxs=0
            if(j==1 || j==Nx || j==Nx1)
                LL(GL)=kqx; CL(GL)=kqx; VL(GL)=1; GL=GL+1;
                R(kqx)=0;
            % Top boundary: free slip/no slip
            elseif(i==1)
                LL(GL)=kqx; CL(GL)=kqx; VL(GL)=1; GL=GL+1;
                LL(GL)=kqx; CL(GL)=kqx+6; VL(GL)=bcftop; GL=GL+1;
                R(kqx)=0;
            else
            % Bottom boundary: free slip/no slip
                LL(GL)=kqx; CL(GL)=kqx-6; VL(GL)=bcfbottom; GL=GL+1;
                LL(GL)=kqx; CL(GL)=kqx; VL(GL)=1; GL=GL+1;
                R(kqx)=0;
            end
        else
        % Internal points: x-Darcy eq.
        % Rx*qxDarcy+dP/dx=RHOfluid*gx
        %     P1-qxD-P2
        % Left part
        LL(GL)=kqx; CL(GL)=kqx; VL(GL)=RX(i,j); GL=GL+1; % qxD
        LL(GL)=kqx; CL(GL)=kpf; VL(GL)=-pscale/dx; GL=GL+1; % P1
        LL(GL)=kqx; CL(GL)=kpf+Ny1*6; VL(GL)=pscale/dx; GL=GL+1; % P2
        % Right part
        R(kqx)=RHOFX(i,j)*gx;
        end
        
        % qyDarcy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Top BC, Bottom BC, Ghost nodes: 1*qys=0
            if(i==1 && j>1 && j<Nx1)
                % %            qyD1
                % %              |
                % %      qxD1--Pt,Pf--qxD2
                % %              |
                % %            qyD2
                % %
                % % Left part
                % LL(GL)=kqy; CL(GL)=kqx+6-Ny1*6; VL(GL)=-1/dx; GL=GL+1; % qxD1
                % LL(GL)=kqy; CL(GL)=kqx+6; VL(GL)=1/dx; GL=GL+1; % qxD2
                % LL(GL)=kqy; CL(GL)=kqy; VL(GL)=-1/dy; GL=GL+1; % qyD1
                % LL(GL)=kqy; CL(GL)=kqy+6; VL(GL)=1/dy; GL=GL+1; % qyD2
                % % Right part
                % R(kqy)=0;
                LL(GL)=kqy; CL(GL)=kqy; VL(GL)=1; GL=GL+1;
                R(kqy)=0;
            elseif(i==Ny || i==Ny1)
                LL(GL)=kqy; CL(GL)=kqy; VL(GL)=1; GL=GL+1;
                R(kqy)=0;
            % Left boundary: free slip/no slip
            elseif(j==1)
                LL(GL)=kqy; CL(GL)=kqy; VL(GL)=1; GL=GL+1;
                LL(GL)=kqy; CL(GL)=kqy+6*Ny1; VL(GL)=bcfleft; GL=GL+1;
                R(kqy)=0;
            else
            % Right boundary: free slip/no slip
                LL(GL)=kqy; CL(GL)=kqy-6*Ny1; VL(GL)=bcfright; GL=GL+1;
                LL(GL)=kqy; CL(GL)=kqy; VL(GL)=1; GL=GL+1;
                R(kqy)=0;
            end
        else
        % Internal points: y-Darcy eq.
        % Rx*qxDarcy+dP/dx=RHOfluid*gx
        %      P1
        %      |
        %     qxD
        %      |
        %      P2
        % Left part
        LL(GL)=kqy; CL(GL)=kqy; VL(GL)=RY(i,j); GL=GL+1; % qyD
        LL(GL)=kqy; CL(GL)=kpf; VL(GL)=-pscale/dy; GL=GL+1; % P1
        LL(GL)=kqy; CL(GL)=kpf+6; VL(GL)=pscale/dy; GL=GL+1; % P2
        % Right part
        R(kqy)=RHOFY(i,j)*gy;
        end
        
        % Pfluid equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            % Ghost nodes: 1*P=0
            LL(GL)=kpf; CL(GL)=kpf; VL(GL)=1*pscale; GL=GL+1;
            R(kpf)=0; % Right part
        % % Real BC
        % elseif(i==2)
        %     LL(GL)=kpf; CL(GL)=kpf; VL(GL)=1*pscale; GL=GL+1;
        %     R(kpf)=ptop; % Right part
        else
        % Internal points: continuity eq.
        % dqxD/dx+dqyD/dy-(Ptotal-Pfluid)/ETHAphi/(1-PHI)-Kbw*BETTAd*(DsPt/dt-1/Ksk*DfPf/dt)=PHI*(ALPHAf*DfT/Dt-ALPHAs*DsT/Dt)
        %            qyD1
        %              |
        %      qxD1--Pt,Pf--qxD2
        %              |
        %            qyD2
        %
        % Pores compressibility
        BETTAPHI=0;
        % Drained compressibility
        BETTADRAINED=(BETTAPHI+BETTASOLID(i,j))/(1-PHI(i,j)); 
        % Biott-Willis koefficient
        KBW=1-BETTASOLID(i,j)/BETTADRAINED; 
        % 1/Skempton_koefficient
        if((BETTADRAINED-BETTASOLID(i,j))>0)
            KSK1=(BETTADRAINED-BETTASOLID(i,j)+PHI(i,j)*(BETTAFLUID(i,j)-BETTASOLID(i,j)))/(BETTADRAINED-BETTASOLID(i,j));
        else
            KSK1=0;
        end
        % Left part
        LL(GL)=kpf; CL(GL)=kqx-Ny1*6; VL(GL)=-1/dx; GL=GL+1; % qxD1
        LL(GL)=kpf; CL(GL)=kqx; VL(GL)=1/dx; GL=GL+1; % qxD2
        LL(GL)=kpf; CL(GL)=kqy-6; VL(GL)=-1/dy; GL=GL+1; % qyD1
        LL(GL)=kpf; CL(GL)=kqy; VL(GL)=1/dy; GL=GL+1; % qyD2
        LL(GL)=kpf; CL(GL)=kpm; VL(GL)=-pscale*(ETAPHI(i,j)/(1-PHI(i,j))+gggbkoef*BETTADRAINED*KBW/dt); GL=GL+1; % qyD2
        LL(GL)=kpf; CL(GL)=kpf; VL(GL)=pscale*(ETAPHI(i,j)/(1-PHI(i,j))+gggbkoef*BETTADRAINED*KBW*KSK1/dt); GL=GL+1; % qyD2
        % Right part
        R(kpf)=gggbkoef*(-BETTADRAINED*KBW*(PT0(i,j)-KSK1*PF0(i,j))/dt+PHI(i,j)*(ALPHAFLUID(i,j)*(tk2(i,j)-tkf1(i,j))/dt-ALPHASOLID(i,j)*(tk2(i,j)-tk1(i,j))/dt));
        end
        
    end
end
%%
% 4) Solving matrixes, reloading solution
% S=L\R; % Obtaining algebraic vector of solutions S()
GLcur=GL-1
% MAKE SPARSE MATRIX
L=sparse(LL(1:GL-1),CL(1:GL-1),VL(1:GL-1));

A=1
% 6) Solving matrix
S=L\R;
% S=mldivide(L,R);
B=1


% Reload solutions S() to vx(), vy(), p()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % P fluid
        % Reload solution
        vx(i,j)=S(kvx);
        vy(i,j)=S(kvy);
        pr(i,j)=S(kpm)*pscale;
        qxD(i,j)=S(kqx);
        qyD(i,j)=S(kqy);
        pf(i,j)=S(kpf)*pscale;
    end
end


dprcorr=(ptop-pr(2,2))
pr=pr+dprcorr;
pf=pf+dprcorr;
pf0=pf0+dprcorr;
PF0=PF0+dprcorr;
PT0=PT0+dprcorr;
aaa(1,1)=vx(5,7);
aaa(2,1)=vy(5,7);
aaa(3,1)=pr(5,7);
aaa(4,1)=qxD(5,7);
aaa(5,1)=qyD(5,7);
aaa(6,1)=pf(5,7);

% Compute Stress and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Strain rate SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
SXX=zeros(Ny1,Nx1); % Strain rate SIGMAxx, Pa
DIVV=zeros(Ny1,Nx1); % div(v)
for j=2:1:Nx
    for i=2:1:Ny
        % DIVV
        DIVV(i,j)=(vx(i,j)-vx(i,j-1))/dx+(vy(i,j)-vy(i-1,j))/dy;
        % EXX
        EXX(i,j)=((vx(i,j)-vx(i,j-1))/dx-(vy(i,j)-vy(i-1,j))/dy)/2;
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j);
    end
end

% Compute Dln((1-PHI)/PHI)/Dt, DPHI/Dt
APHI=zeros(Ny1,Nx1);
APHI0=zeros(Ny1,Nx1);
aphimax=0;
for j=2:1:Nx
    for i=2:1:Ny
        % Compute Dln((1-PHI)/PHI)/Dt
        if(PHI(i,j)>0)
            APHI(i,j)=((pr(i,j)-pf(i,j))/(1-PHI(i,j))*ETAPHI(i,j))/PHI(i,j);
            aphimax=max(aphimax,abs(APHI(i,j)));
        end
        % Compute DPHI/Dt
        APHI0(i,j)=-(pr(i,j)-pf(i,j))*ETAPHI(i,j);
    end
end

% Apply Symmetry
% Top
APHI(1,2:Nx)=APHI(2,2:Nx);    
APHI0(1,2:Nx)=APHI0(2,2:Nx);    
PHI(1,2:Nx)=PHI(2,2:Nx);    
pr(1,2:Nx)=pr(2,2:Nx);    
pf(1,2:Nx)=pf(2,2:Nx);    
EXX(1,2:Nx)=EXX(2,2:Nx);    
SXX(1,2:Nx)=SXX(2,2:Nx);    
% Bottom
APHI(Ny1,2:Nx)=APHI(Ny,2:Nx);    
APHI0(Ny1,2:Nx)=APHI0(Ny,2:Nx);    
PHI(Ny1,2:Nx)=PHI(Ny,2:Nx);    
pr(Ny1,2:Nx)=pr(Ny,2:Nx);    
pf(Ny1,2:Nx)=pf(Ny,2:Nx);    
EXX(Ny1,2:Nx)=EXX(Ny,2:Nx);    
SXX(Ny1,2:Nx)=SXX(Ny,2:Nx);    
% Left
APHI(:,1)=APHI(:,2);    
APHI0(:,1)=APHI0(:,2);    
PHI(:,1)=PHI(:,2);    
pr(:,1)=pr(:,2);    
pf(:,1)=pf(:,2);    
EXX(:,1)=EXX(:,2);    
SXX(:,1)=SXX(:,2);    
% Right
APHI(:,Nx1)=APHI(:,Nx);    
APHI0(:,Nx1)=APHI0(:,Nx);    
PHI(:,Nx1)=PHI(:,Nx);    
pr(:,Nx1)=pr(:,Nx);    
pf(:,Nx1)=pf(:,Nx); 
EXX(:,Nx1)=EXX(:,Nx); 
SXX(:,Nx1)=SXX(:,Nx); 

% Compute solid pressure
ps=(pr-pf.*PHI)./(1-PHI);

% Compute fluid velocity
% Vx fluid
for j=1:1:Nx
    for i=2:1:Ny
        vxf(i,j)=0;
        if(PHIX(i,j)>0)
            vxf(i,j)=qxD(i,j)/PHIX(i,j);
        end
    end
end
% Apply BC
% Top
vxf(1,:)=-bcftop*vxf(2,:);    
% Bottom
vxf(Ny1,:)=-bcfbottom*vxf(Ny,:);    
% Vy fluid
for j=2:1:Nx
    for i=1:1:Ny
        vyf(i,j)=0;
        if(PHIY(i,j)>0)
            vyf(i,j)=qyD(i,j)/PHIY(i,j);
        end
    end
end
% Apply BC
% Left
vyf(:,1)=-bcfleft*vyf(:,2);    
% Right
vyf(:,Nx1)=-bcfright*vyf(:,Nx);     
% Add solid velocity
vxf0=vxf; vxf=vxf+vx;
vyf0=vyf; vyf=vyf+vy;

% Compute shear heating HS in pressure nodes
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)*EXY(i,j)+SXY(i-1,j)*EXY(i-1,j)+...
            SXY(i,j-1)*EXY(i,j-1)+SXY(i-1,j-1)*EXY(i-1,j-1))/4;
        % HS
        HS(i,j)=2*SXX(i,j)*EXX(i,j)+2*SXYEXY+...
            (pr(i,j)-pf(i,j))^2/(1-PHI(i,j))*ETAPHI(i,j)+...
            (RX(i,j-1)*qxD(i,j-1)^2+RX(i,j)*qxD(i,j)^2)/2+...
            (RY(i-1,j)*qyD(i-1,j)^2+RY(i,j)*qyD(i,j)^2)/2;
    end
end
%%
% Compute total adiabatic heating HA*dt in pressure nodes
HA=zeros(Ny1,Nx1); % Shear heating, W/m^3
if(timestep>1)
for j=2:1:Nx
    for i=2:1:Ny
        % HA
        % % Average vy, vx, vxf, vyf
        % VXP=(vx(i,j)+vx(i,j-1))/2;
        % VYP=(vy(i,j)+vy(i-1,j))/2;
        % VXFP=(vxf(i,j)+vxf(i,j-1))/2;
        % VYFP=(vyf(i,j)+vyf(i-1,j))/2;
        % % Evaluate DPsolid/Dt with upwind differences
        % if(VXP<0)
        %     dpsdx=(ps(i,j)-ps(i,j-1))/dx;
        %     dprdx=(pr(i,j)-pr(i,j-1))/dx;
        % else
        %     dpsdx=(ps(i,j+1)-ps(i,j))/dx;
        %     dprdx=(pr(i,j+1)-pr(i,j))/dx;
        % end
        % if(VYP<0)
        %     dpsdy=(ps(i,j)-ps(i-1,j))/dy;
        %     dprdy=(pr(i,j)-pr(i-1,j))/dy;
        % else
        %     dpsdy=(ps(i+1,j)-ps(i,j))/dy;
        %     dprdy=(pr(i+1,j)-pr(i,j))/dy;
        % end
        % dpsdt=VXP*dpsdx+VYP*dpsdy;
        % dprdt=VXP*dprdx+VYP*dprdy;
        % % Evaluate DPfluid/Dt with upwind differences
        % if(VXFP>0)
        %     dpfdx=(pf(i,j)-pf(i,j-1))/dx;
        % else
        %     dpfdx=(pf(i,j+1)-pf(i,j))/dx;
        % end
        % if(VYFP>0)
        %     dpfdy=(pf(i,j)-pf(i-1,j))/dy;
        % else
        %     dpfdy=(pf(i+1,j)-pf(i,j))/dy;
        % end
        % dpfdt=VXFP*dpsdx+VYFP*dpsdy;
        % HA
        % HA(i,j)=(1-PHI(i,j))*tk1(i,j)*ALPHA(i,j)*dpsdt+...
        %     PHI(i,j)*tk1(i,j)*alphafluid*dpfdt;
        % HA(i,j)=ALPHA(i,j)*dprdt;
        % HA(i,j)=ALPHA(i,j)*((1-PHI(i,j))*dpsdt+PHI(i,j)*dpfdt);
        dprdt=(pr(i,j)-PT00(i,j))/dt;
        HA(i,j)=ALPHA(i,j)*dprdt;
    end
end
end


%%

% Define timestep
if(titer==1)
    dt=dt*dtkoef;
    % Solid velocity
    maxvx=max(max(abs(vx)));
    maxvy=max(max(abs(vy)));
    if(dt*maxvx>dxymax*dx)
        dt=dxymax*dx/maxvx;
    end
    if(dt*maxvy>dxymax*dy)
        dt=dxymax*dy/maxvy;
    end
    % Fluid velocity
    maxvxf=max(max(abs(vxf)));
    maxvyf=max(max(abs(vyf)));
    if(dt*maxvxf>dxymax*dx)
        dt=dxymax*dx/maxvxf;
    end
    if(dt*maxvyf>dxymax*dy)
        dt=dxymax*dy/maxvyf;
    end
    % Porosity change
    if(aphimax*dt>dphimax)
        dt=dphimax/aphimax;
    end
    dt
end

%%
% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
GLT=1;
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % External points
        if(i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % Top BC: T=tktop
            if(i==1 && j>1 && j<Nx1)
                LLT(GLT)=gk; CLT(GLT)=gk; VLT(GLT)=1; GLT=GLT+1;
                LLT(GLT)=gk; CLT(GLT)=gk+1; VLT(GLT)=1; GLT=GLT+1;
                RT(gk)=tktop*2; % Right part
            end
            % Bottom BC: T=1500
            if(i==Ny1 && j>1 && j<Nx1)
                LLT(GLT)=gk; CLT(GLT)=gk-1; VLT(GLT)=1; GLT=GLT+1;
                LLT(GLT)=gk; CLT(GLT)=gk; VLT(GLT)=1; GLT=GLT+1;
                RT(gk)=tkbottom*2; % Right part
            end
            % Left BC: dT/dx=0
            if(j==1)
                LLT(GLT)=gk; CLT(GLT)=gk; VLT(GLT)=1; GLT=GLT+1;
                LLT(GLT)=gk; CLT(GLT)=gk+Ny1; VLT(GLT)=-1; GLT=GLT+1;
                RT(gk)=0; % Right part
            end
            % Right BC: dT/dx=0
            if(j==Nx1)
                LLT(GLT)=gk; CLT(GLT)=gk-Ny1; VLT(GLT)=-1; GLT=GLT+1;
                LLT(GLT)=gk; CLT(GLT)=gk; VLT(GLT)=1; GLT=GLT+1;
                RT(gk)=0; % Right part
            end
        else
        % Internal points: Temperature eq.
        % RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
        %          Tdt2
        %           |
        %          Ky1
        %           |
        %Tdt1-Kx1-T03,Tdt3-Kx2-Tdt5
        %           |
        %          Ky2
        %           |
        %          Tdt4
        %
        % Left part
        Kx1=KX(i,j-1); 
        Kx2=KX(i,j); 
        Ky1=KY(i-1,j); 
        Ky2=KY(i,j); 
        LLT(GLT)=gk; CLT(GLT)=gk-Ny1; VLT(GLT)=-Kx1/dx^2; GLT=GLT+1; % T1
        LLT(GLT)=gk; CLT(GLT)=gk-1; VLT(GLT)=-Ky1/dy^2; GLT=GLT+1; % T2
        LLT(GLT)=gk; CLT(GLT)=gk; VLT(GLT)=RHOCP(i,j)/dt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; GLT=GLT+1; % T3
        LLT(GLT)=gk; CLT(GLT)=gk+1; VLT(GLT)=-Ky2/dy^2; GLT=GLT+1; % T4
        LLT(GLT)=gk; CLT(GLT)=gk+Ny1; VLT(GLT)=-Kx2/dx^2; GLT=GLT+1; % T5
        % Right part
        RT(gk)=RHOCP(i,j)/dt*tk1(i,j)+HR(i,j)+HA(i,j)+HS(i,j)+DHP(i,j);
        end
    end
end
%%
% Solving matrixes
% ST=LT\RT; % Obtaining algebraic vector of solutions ST()
GLTcur=GLT-1
% MAKE SPARSE MATRIX
LT=sparse(LLT(1:GLT-1),CLT(1:GLT-1),VLT(1:GLT-1));

% 6) Solving matrix
ST=LT\RT;
% ST=mldivide(LT,RT);


% Reload solutions ST() to geometrical array Tdt()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        tk2(i,j)=ST(gk);
    end
end
aaa(7,1)=tk2(24,20);
%%
% Compute DT
DT=tk2-tk1;
titer;
dt;

% pause(100)

% Apply thermal timestepping condition
if(titer==1)
    maxDTcurrent=max(max(abs(DT)));
    if(maxDTcurrent>DTmax)
        dt=dt/maxDTcurrent*DTmax;
    end
    dt;
end

if(1==0)
figure(3);colormap('Jet');clf
subplot(4,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['logETA,Pa*s timestep=',num2str(timestep)])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title(['Ptotal,Pa time(My)=',num2str(timesum/1e+6/(365.25*24*3600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,3)
pcolor(xvx/1000,yvx/1000,vx)
shading interp;
axis ij image;
colorbar
title('vx, m/s')
title(['vx,m/s dt(y)=',num2str(dt/(365.25*24*3600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,4)
pcolor(xvy/1000,yvy/1000,vy)
shading interp;
axis ij image;
colorbar
title('vy, m/s')
title(['vy, m/s iter=',num2str(titer)])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,5)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('HS, W/m^3')

subplot(4,4,6)
pcolor(xp/1000,yp/1000,RHO)
shading interp;
axis ij image;
colorbar
title('RHO, kg/m^3')

subplot(4,4,7)
pcolor(xp/1000,yp/1000,pf-pf0)
shading interp;
axis ij image;
colorbar
title('pf-pf0,Pa')

subplot(4,4,8)
pcolor(xvx/1000,yvx/1000,log10(KX))
shading interp;
axis ij image;
colorbar
title('logK, W/m/K')

subplot(4,4,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')

subplot(4,4,10)
pcolor(xp/1000,yp/1000,PHI)
shading interp;
axis ij image;
colorbar
title('PHI')

subplot(4,4,11)
pcolor(xp/1000,yp/1000,pr-pf)
shading interp;
axis ij image;
colorbar
title('Ptotal-Pfluid, Pa')

subplot(4,4,12)
pcolor(xvx/1000,yvx/1000,qxD)
shading interp;
axis ij image;
colorbar
title('qxDarcy, m/s')

subplot(4,4,13)
pcolor(xvy/1000,yvy/1000,qyD)
shading interp;
axis ij image;
colorbar
title('qyDarcy, m/s')


subplot(4,4,14)
pcolor(xvx/1000,yvx/1000,XHT)
shading interp;
axis ij image;
colorbar
title('XHT')


subplot(4,4,15)
pcolor(xp/1000,yp/1000,log10(1./(ETAPHI+1e-30)))
shading interp;
axis ij image;
colorbar
title('logETAphi, Pa*s')

% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,APHI)
% shading interp;
% axis ij image;
% colorbar
% title('Dln((1-PHI)/PHI)/dt')
subplot(4,4,16)
pcolor(xp/1000,yp/1000,DMP)
shading interp;
axis ij image;
colorbar
title('DMP, 1/s')
% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,DHP)
% shading interp;
% axis ij image;
% colorbar
% title('DHP, J/m^3/s')


% figure(2);colormap('Jet');clf
% plot(xm,ym,'. k');
% axis ij image


figure(5);colormap('Jet');clf
subplot(2,3,1)
pcolor(xvx/1000,yvx/1000,log10(RX))
shading interp;
axis ij image;
colorbar
title('logRX')

subplot(2,3,2)
pcolor(xvy/1000,yvy/1000,log10(RY))
shading interp;
axis ij image;
colorbar
title('logRY')

subplot(2,3,3)
pcolor(xvx/1000,yvx/1000,RHOFX)
shading interp;
axis ij image;
colorbar
title('FROFX')

subplot(2,3,4)
pcolor(xvy/1000,yvy/1000,RHOFY)
shading interp;
axis ij image;
colorbar
title('RHOFY')

subplot(2,3,5)
pcolor(xp/1000,yp/1000,log10(1./(ETAPHI+1e-30)))
shading interp;
axis ij image;
colorbar
title('logETAPHI')

pause(0.01)
end


%
dt
% Stop iterations
pferrcur=max(max(abs(pf-pf0)))
DMPmax=max(max(abs(DMP)));
if((pferrcur<pferrmax && (titer>2 || DMPmax<=0)) || (titer>=titermax1 && timestep<=titerstep1))
    break;
end

end
DT0=DT;


% Apply subgrid diffusion on markers
if(dsubgridt>0)
TKSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node T difference
    dtkm0=tkm(m)-(tk1(i,j)*wtmij+tk1(i+1,j)*wtmi1j+...
            tk1(i,j+1)*wtmij1+tk1(i+1,j+1)*wtmi1j1);
    % Compute marker parameters
    if(tm(m)<4)
        % Rocks
         ktotalm=(ksolidm(tm(m))*kfluid/2+((ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))^2)/16)^0.5-(ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))/4;
    else
        % Sticky air
        ktotalm=ksolidm(tm(m));
    end    % Relax temperature difference
    dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dt/rhocptotalm(m)*(2/dx^2+2/dy^2));
    % Correct marker temperature
    ddtkm=dtkm1-dtkm0;
    tkm(m)=tkm(m)+ddtkm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TKSUM(i,j)=TKSUM(i,j)+ddtkm*rhocptotalm(m)*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocptotalm(m)*wtmij;
    % i+1,j Node
    TKSUM(i+1,j)=TKSUM(i+1,j)+ddtkm*rhocptotalm(m)*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocptotalm(m)*wtmi1j;
    % i,j+1 Node
    TKSUM(i,j+1)=TKSUM(i,j+1)+ddtkm*rhocptotalm(m)*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocptotalm(m)*wtmij1;
    % i+1,j+1 Node
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+ddtkm*rhocptotalm(m)*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocptotalm(m)*wtmi1j1;
end
end
% Compute DTsubgrid
DTsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOCPSUM(i,j)>0)
            DTsubgrid(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Correct DT
DT=DT-DTsubgrid;
end

% Interpolate DT to markers
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    tkm(m)=tkm(m)+DT(i,j)*wtmij+DT(i+1,j)*wtmi1j+...
            DT(i,j+1)*wtmij1+DT(i+1,j+1)*wtmi1j1;
    % Interpolate tk2 at 1st timestep
    if(timestep==1)
        tkm(m)=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    end
end
end


% Update porosity for (de)compaction
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize && tm(m)<4)
    phim(m)=min(max(phicomm(m),phimin),phimax);
    % Interpolate APHI
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute Dln((1-PHI)/PHI)/Dt
    aphim=APHI(i,j)*wtmij+APHI(i+1,j)*wtmi1j+...
        APHI(i,j+1)*wtmij1+APHI(i+1,j+1)*wtmi1j1;
    aphim0=APHI0(i,j)*wtmij+APHI0(i+1,j)*wtmi1j+...
        APHI0(i,j+1)*wtmij1+APHI0(i+1,j+1)*wtmi1j1;
    phimm=PHI(i,j)*wtmij+PHI(i+1,j)*wtmi1j+...
        PHI(i,j+1)*wtmij1+PHI(i+1,j+1)*wtmi1j1;
    % % Change Porosity

    phikf=1/((1-phim(m))*exp(aphim*dt)+phim(m));

    dphim=phim(m)*(phikf-1);

    phicomm(m)=min(max(phicomm(m)+dphim,0),1);
    phim(m)=min(max(phicomm(m),phimin),phimax);
end
end

%UPDATE HYDRATION
for j = 2:Nx
    for i = 2:Ny
        DIVH(i,j) = ( RHOWATER(i,j)   * qxD(i,j)   - RHOWATER(i,j-1) * qxD(i,j-1) ) / dx ...
                   + (RHOWATER(i,j)   * qyD(i,j)   - RHOWATER(i-1,j) * qyD(i-1,j) ) / dy;
    end
end
%B.C
DIVH(1,2:Nx)   = DIVH(2,2:Nx);
DIVH(Ny1,2:Nx) = DIVH(Ny,2:Nx);
DIVH(:,1)      = DIVH(:,2);
DIVH(:,Nx1)    = DIVH(:,Nx);
for m = 1:1:marknum
if (xm(m) >= 0 && xm(m) <= xsize && ym(m) >= 0 && ym(m) <= ysize && tm(m) < 4)
    % Compute marker parameters
    % Compute hydration of solid
    % Interpolate fluid pressure and temperature to the marker
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;

    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Interpolate nodal total pressure and temperature to the marker
    ptnm=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
    tknm=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    divH=DIVH(i,j)*wtmij+DIVH(i+1,j)*wtmi1j+...
            DIVH(i,j+1)*wtmij1+DIVH(i+1,j+1)*wtmi1j1;
    rho_water=RHOWATER(i,j)*wtmij+RHOWATER(i+1,j)*wtmi1j+...
            RHOWATER(i,j+1)*wtmij1+RHOWATER(i+1,j+1)*wtmi1j1;
    rho_solid=RHOSOLID(i,j)*wtmij+RHOSOLID(i+1,j)*wtmi1j+...
            RHOSOLID(i,j+1)*wtmij1+RHOSOLID(i+1,j+1)*wtmi1j1;

    n_tot = (1 - phicomm(m)) * rho_solid / Mdry + phicomm(m) * rho_water / Mf;

    dXH_darcy = - divH / (n_tot * Mf);

    XHtotalm(m) = XHtotalm(m) + dXH_darcy * dt;
    XHtotalm(m) = max(0, min(1, XHtotalm(m)));

    % Thermodynamic computations
    T=tknm; % Temperature, K
    if(T<=0)
        T=tkm(m);
    end
    PB=(ptnm)*1e-9; % Pressure, Pa=1e+9 GPa
    Wth20=XHtotalm(m);

    nt=fix((T-tkmin)/tkstp)+1;
    np=fix((PB-pbmin)/pbstp)+1;
    nwt=fix((Wth20-wth20min)/Wth20stp)+1;

    if(nt<1) 
        nt=1; 
    elseif(nt>tknum-1) 
        nt=tknum-1;
    end
    if(np<1) 
        np=1;
    elseif(np>pbnum-1) 
        np=pbnum-1;
    end
    if(nwt<1) 
        nwt=1;
    elseif(nwt>nwtnum-1) 
        nwt=nwtnum-1;
    end
    % Compute relative distances from the marker to the P-T-node
    dtt=(T-tkj(nt))/tkstp;
    dp=(PB-pbi(np))/pbstp;
    dwt=(Wth20-wth20k(nwt))/Wth20stp;

    w000 = (1-dp)*(1-dtt)*(1-dwt);
    w100 =    dp *(1-dtt)*(1-dwt);
    w010 = (1-dp)*   dtt *(1-dwt);
    w110 =    dp *   dtt *(1-dwt);
    w001 = (1-dp)*(1-dtt)*   dwt;
    w101 =    dp *(1-dtt)*   dwt;
    w011 = (1-dp)*   dtt *   dwt;
    w111 =    dp *   dtt *   dwt;
    % Corner concentrations
    c000 = XH2Osolid(np,  nt,   nwt  );
    c100 = XH2Osolid(np+1,nt,   nwt  );
    c010 = XH2Osolid(np,  nt+1, nwt  );
    c110 = XH2Osolid(np+1,nt+1, nwt  );
    c001 = XH2Osolid(np,  nt,   nwt+1);
    c101 = XH2Osolid(np+1,nt,   nwt+1);
    c011 = XH2Osolid(np,  nt+1, nwt+1);
    c111 = XH2Osolid(np+1,nt+1, nwt+1);

    XHsolidm0(m) = w000*c000 + w100*c100 + w010*c010 + w110*c110 + ...
                   w001*c001 + w101*c101 + w011*c011 + w111*c111;
    XBsolidm0(m)=(1-XHsolidm0(m));
end
end

% Compute viscosity on markers
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize && tm(m)<4)
    % Interpolate EXY for the marker location
    % Interpolation from basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute EXY
    exym=EXY(i,j)*wtmij+EXY(i+1,j)*wtmi1j+...
            EXY(i,j+1)*wtmij1+EXY(i+1,j+1)*wtmi1j1;
    
    % Interpolate P,T,EXX for the marker location
    % Interpolation from P nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute P,T,EXX
    P=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
    T=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    exxm=EXX(i,j)*wtmij+EXX(i+1,j)*wtmi1j+...
            EXX(i,j+1)*wtmij1+EXX(i+1,j+1)*wtmi1j1;
    % Compute second strain rate invariant
    epsin=(exxm^2+exym^2)^0.5;
    % Compute marker viscosity 
    % Diffusion creep
    eta_dif=etamax;

    % Dislocation creep
    % Compute and check activation energy exponent
    eterm=(Ea_dis(flow_law)+P*Va_dis(flow_law))/8.314/T ;
    if (eterm>150) 
        eterm=150; 
    end
    eterm=exp(-eterm);

    % Compute initial approximation for stress invariant
    etacur0=etamin;
    sigin0=2*epsin*etacur0;
    % Compute dislocation creep viscosity
    eta_dis=sigin0/(2*A_dis(flow_law)*(sigin0*1e-6)^n_dis(flow_law)*eterm);
    % Compute effective viscosity for dry creep
    etadry1=1/(1/eta_dis+1/eta_dif);
    % Check if to continue
    if(etadry1<=etamin)
        etasolidm(m)=etamin;
    else
        % Compute second approximation for stress invariant
        etacur1=etamax;
        sigin1=2*epsin*etacur1;
        % Compute dislocation creep viscosity
        eta_dis=sigin1/(2*A_dis(flow_law)*(sigin1*1e-6)^n_dis(flow_law)*eterm);
        % Compute effective viscosity for dry creep
        etadry1=1/(1/eta_dis+1/eta_dif);
        % Check if to continue
        if(etadry1>=etamax)
            etasolidm(m)=etamax;
        else
            % Bisection
            for nnn=1:1:100
                % Compute current approximation for stress invariant
                etacur2=exp((log(etacur0*etacur1))/2);
                sigin2=2*epsin*etacur2;
                % Compute dislocation creep viscosity
                eta_dis=sigin2/(2*A_dis(flow_law)*(sigin2*1e-6)^n_dis(flow_law)*eterm);
                % Compute effective viscosity for dry creep
                etadry1=1/(1/eta_dis+1/eta_dif);
                % Exit
                if(abs(log(etadry1/etacur2))<1e-7)
                    etasolidm(m)=etacur2;
                    break;
                end
                % Bisection
                if(etadry1>etacur2)
                    etacur0=etacur2;
                else
                    etacur1=etacur2;
                end
            end
        end
    end
end
end



% Compute solid velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxp(i,j)=(vx(i,j)+vx(i,j-1))/2;
    end
end
% Apply BC
% Top
vxp(1,2:Nx-1)=-bctop*vxp(2,2:Nx-1);    
% Bottom
vxp(Ny1,2:Nx-1)=-bcbottom*vxp(Ny,2:Nx-1);    
% Left
vxp(:,1)=-vxp(:,2);
% Right
vxp(:,Nx1)=-vxp(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vyp(i,j)=(vy(i,j)+vy(i-1,j))/2;
    end
end    
% Apply BC
% Left
vyp(2:Ny-1,1)=-bcleft*vyp(2:Ny-1,2);    
% Right
vyp(2:Ny-1,Nx1)=-bcright*vyp(2:Ny-1,Nx);     
% Top
vyp(1,:)=-vyp(2,:);
% Bottom
vyp(Ny1,:)=-vyp(Ny,:);

% Compute fluid velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxpf(i,j)=(vxf0(i,j)+vxf0(i,j-1))/2;
    end
end
% Apply BC
% Top
vxpf(1,2:Nx-1)=-bcftop*vxpf(2,2:Nx-1);    
% Bottom
vxpf(Ny1,2:Nx-1)=-bcfbottom*vxpf(Ny,2:Nx-1);    
% Left
vxpf(:,1)=-vxpf(:,2);
% Right
vxpf(:,Nx1)=-vxpf(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vypf(i,j)=(vyf0(i,j)+vyf0(i-1,j))/2;
    end
end    
% Apply BC
% Left
vypf(2:Ny-1,1)=-bcfleft*vypf(2:Ny-1,2);    
% Right
vypf(2:Ny-1,Nx1)=-bcfright*vypf(2:Ny-1,Nx);    
% Top
vypf(1,:)=-vypf(2,:);
% Bottom
vypf(Ny1,:)=-vypf(Ny,:);
% Add solid velocity
vxpf=vxpf+vxp;
vypf=vypf+vyp;

% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
    % Save initial marker coordinates
    xA=xm(m);
    yA=ym(m);
    for rk=1:1:4
        % Interpolate vx
        j=fix((xm(m)-xvx(1))/dx)+1;
        i=fix((ym(m)-yvx(1))/dy)+1;
        % Normalized marker distances
        dxmj=(xm(m)-xvx(j));
        dymi=(ym(m)-yvx(i));
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute marker vx
        vxm(rk)=vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1;
        % Compute correction
        dvxmcorr=0;
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                dvxmcorr=1/2*((dxmj/dx-0.5)^2)*...
                     ((1-dymi/dy)*(vx(i,j)-2*vx(i,j+1)+vx(i,j+2))+...
                     (dymi/dy)*(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2)));
            end
        else
            if(j>1)
                dvxmcorr=1/2*((dxmj/dx-0.5)^2)*...
                    ((1-dymi/dy)*(vx(i,j-1)-2*vx(i,j)+vx(i,j+1))+...
                     (dymi/dy)*(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1)));
            end
        end
        % dvxmcorr=0;
        vxm(rk)=vxm(rk)+dvxmcorr;
        
        % Interpolate vy
        j=fix((xm(m)-xvy(1))/dx)+1;
        i=fix((ym(m)-yvy(1))/dy)+1;
        % Normalized marker distances
        dxmj=(xm(m)-xvy(j));
        dymi=(ym(m)-yvy(i));
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute marker vy
        vym(rk)=vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1;
        dvymcorr=0;
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                dvymcorr=1/2*((dymi/dy-0.5)^2)*...
                    ((1-dxmj/dx)*(vy(i,j)-2*vy(i+1,j)+vy(i+2,j))+...
                     (dxmj/dx)*(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1)));
            end      
        else
            if(i>1)
                dvymcorr=1/2*((dymi/dy-0.5)^2)*...
                    ((1-dxmj/dx)*(vy(i-1,j)-2*vy(i,j)+vy(i+1,j))+...
                    (dxmj/dx)*(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1)));
            end
        end
        % dvymcorr=0;
        vym(rk)=vym(rk)+dvymcorr;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xm(m)=xA+dt/2*vxm(rk);
            ym(m)=yA+dt/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dt*vxm(rk);
            ym(m)=yA+dt*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dt*vxmeff;
    ym(m)=ym(m)+dt*vymeff;
end 
end


% Backtracing total and fluid pressure
% Internal P-nodes
for j1=2:1:Nx
for i1=2:1:Ny
    % Backtracing node with solid velocity
    xcur=xp(j1);
    ycur=yp(i1);
    xA=xcur;
    yA=ycur;
    for rk=1:1:4
        % Interpolate vx
        j=fix((xcur-xvx(1))/dx)+1;
        i=fix((ycur-yvx(1))/dy)+1;
        % Normalized marker distances
        dxmj=(xcur-xvx(j));
        dymi=(ycur-yvx(i));
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute marker vx
        vxm(rk)=vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1;
        % Compute correction
        dvxmcorr=0;
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                dvxmcorr=1/2*((dxmj/dx-0.5)^2)*...
                     ((1-dymi/dy)*(vx(i,j)-2*vx(i,j+1)+vx(i,j+2))+...
                     (dymi/dy)*(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2)));
            end
        else
            if(j>1)
                dvxmcorr=1/2*((dxmj/dx-0.5)^2)*...
                    ((1-dymi/dy)*(vx(i,j-1)-2*vx(i,j)+vx(i,j+1))+...
                     (dymi/dy)*(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1)));
            end
        end
        % dvxmcorr=0;
        vxm(rk)=vxm(rk)+dvxmcorr;
        
        % Interpolate vy
        j=fix((xcur-xvy(1))/dx)+1;
        i=fix((ycur-yvy(1))/dy)+1;
        % Normalized marker distances
        dxmj=(xcur-xvy(j));
        dymi=(ycur-yvy(i));
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute marker vy
        vym(rk)=vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1;
        dvymcorr=0;
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                dvymcorr=1/2*((dymi/dy-0.5)^2)*...
                    ((1-dxmj/dx)*(vy(i,j)-2*vy(i+1,j)+vy(i+2,j))+...
                     (dxmj/dx)*(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1)));
            end      
        else
            if(i>1)
                dvymcorr=1/2*((dymi/dy-0.5)^2)*...
                    ((1-dxmj/dx)*(vy(i-1,j)-2*vy(i,j)+vy(i+1,j))+...
                    (dxmj/dx)*(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1)));
            end
        end
        % dvymcorr=0;
        vym(rk)=vym(rk)+dvymcorr;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xcur=xA-dt/2*vxm(rk);
            ycur=yA-dt/2*vym(rk);
        elseif(rk==3)
            xcur=xA-dt*vxm(rk);
            ycur=yA-dt*vym(rk);
        end
    end
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Trace the node backward
    xcur=xA-dt*vxmeff;
    ycur=yA-dt*vymeff;
    % Interpolate total pressure
    % Define i,j indexes for the upper left node
    j=fix((xcur-xp(1))/dx)+1;
    i=fix((ycur-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xcur-xp(j);
    dymi=ycur-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal old total pressure
    PT0(i1,j1)=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
end
end


% Apply Symmetry
% Top
PT0(1,2:Nx)=PT0(2,2:Nx);    
PF0(1,2:Nx)=PF0(2,2:Nx);    
% Bottom
PT0(Ny1,2:Nx)=PT0(Ny,2:Nx);    
PF0(Ny1,2:Nx)=PF0(Ny,2:Nx);    
% Left
PT0(:,1)=PT0(:,2);    
PF0(:,1)=PF0(:,2);    
% Right
PT0(:,Nx1)=PT0(:,Nx);    
PF0(:,Nx1)=PF0(:,Nx); 


% Add markers to empty areas
marknumold=marknum;
mdis=1e30*ones(Nym,Nxm);
mnum=zeros(Nym,Nxm);
mtyp=zeros(Nym,Nxm);
mpor=zeros(Nym,Nxm);
xxm=dxm/2:dxm:xsize-dxm/2;
yym=dym/2:dym:ysize-dym/2;
for m=1:1:marknum
if(xm(m)>=0 && xm(m)<=xsize && ym(m)>=0 && ym(m)<=ysize)
    
    % Check markers with the nearest nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xxm(1))/dxm)+1;
    i=fix((ym(m)-yym(1))/dym)+1;
    if(j<1)
        j=1;
    elseif(j>Nxm-1)
        j=Nxm-1;
    end
    if(i<1)
        i=1;
    elseif(i>Nym-1)
        i=Nym-1;
    end
    
    % Check nodes
    % i,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i);
    dismij=(dxmj^2+dymi^2)^0.5;
    if(dismij<mdis(i,j))
        mdis(i,j)=dismij;
        mnum(i,j)=m;
        mtyp(i,j)=tm(m);
        mpor(i,j)=phicomm(m);
    end
    % i+1,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i+1);
    dismi1j=(dxmj^2+dymi^2)^0.5;
    if(dismi1j<mdis(i+1,j))
        mdis(i+1,j)=dismi1j;
        mnum(i+1,j)=m;
        mtyp(i+1,j)=tm(m);
        mpor(i+1,j)=phicomm(m);
    end
    % i,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i);
    dismij1=(dxmj^2+dymi^2)^0.5;
    if(dismij1<mdis(i,j+1))
        mdis(i,j+1)=dismij1;
        mnum(i,j+1)=m;
        mtyp(i,j+1)=tm(m);
        mpor(i,j+1)=phicomm(m);
    end
    % i+1,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i+1);
    dismi1j1=(dxmj^2+dymi^2)^0.5;
    if(dismi1j1<mdis(i+1,j+1))
        mdis(i+1,j+1)=dismi1j1;
        mnum(i+1,j+1)=m;
        mtyp(i+1,j+1)=tm(m);
        mpor(i+1,j+1)=phicomm(m);
    end
end
end

dii=5*Nxmc;
djj=5*Nymc;

for j=1:1:Nxm
    for i=1:1:Nym
        if(mnum(i,j)==0)
            % Serch surrounding nodes
            for jj=j-djj:1:j+djj
                for ii=i-dii:1:i+dii
                    if(ii>=1 && ii<=Nym && jj>=1 && jj<=Nxm && mnum(ii,jj)>0)
                        % Compute distance
                        m=mnum(ii,jj);
                        dxmj=xm(m)-xxm(j);
                        dymi=ym(m)-yym(i);
                        dismij=(dxmj^2+dymi^2)^0.5;
                        if(dismij<mdis(i,j))
                            mdis(i,j)=dismij;
                            mnum(i,j)=-m;
                            mtyp(i,j)=tm(m);
                            mpor(i,j)=phicomm(m);
                        end
                    end
                end
            end
            % Add New marker
            if(mnum(i,j)<0)
                % Add marker number
                marknum=marknum+1;
                % Assign marker coordinates
                xm(marknum)=xxm(j)+(rand-0.5)*dxm;
                ym(marknum)=yym(i)+(rand-0.5)*dym;
                % Copy marker properties
                m=-mnum(i,j);
                tm(marknum)=tm(m); % Material type
                tkm(marknum)=tkm(m); % Temperature
                phim(marknum)=phim(m);
                phicomm(marknum)=phicomm(m);
                phinewm(marknum)=phinewm(m);
                etasolidm(marknum)=etasolidm(m);
                dtkfsm(marknum)=dtkfsm(m);
                rhocptotalm(marknum)=rhocptotalm(m);
                XHsolidm(marknum)=XHsolidm(m);
                XBsolidm0(marknum)=XBsolidm0(m);
                XHsolidm0(marknum)=XHsolidm0(m);
                XHtotalm(marknum)=XHtotalm(m);
                pfm0(marknum)=pfm0(m);
                dxbf0(marknum)=dxbf0(m);
                dxhf0(marknum)=dxhf0(m);
            end
        end
    end
end           
marknumnew=marknum;


% Update timesum
timesum=timesum+dt;


if(desktopyn==1)
figure(1);colormap('Jet');clf
subplot(4,5,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['logETA,Pa*s timestep=',num2str(timestep)])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,5,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title(['Ptotal,Pa time(My)=',num2str(timesum/1e+6/(365.25*24*3600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,5,3)
pcolor(xvx/1000,yvx/1000,vx)
shading interp;
axis ij image;
colorbar
title('vx, m/s')
title(['vx,m/s dt(y)=',num2str(dt/(365.25*24*3600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,5,4)
pcolor(xvy/1000,yvy/1000,vy)
shading interp;
axis ij image;
colorbar
title('vy, m/s')
title(['vy, m/s iter=',num2str(titer)])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,5,5)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('HS, W/m^3')

subplot(4,5,6)
pcolor(xp/1000,yp/1000,log10(1./(ETAPHI+1e-30)))
shading interp;
axis ij image;
colorbar
title('logETAphi, Pa*s')

subplot(4,5,7)
pcolor(xp/1000,yp/1000,RHO)
shading interp;
axis ij image;
colorbar
title('RHOtotal, kg/m^3')

subplot(4,5,8)
pcolor(xvx/1000,yvx/1000,log10(KX))
shading interp;
axis ij image;
colorbar
title('logK, W/m/K')

subplot(4,5,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')

subplot(4,5,10)
pcolor(xp/1000,yp/1000,PHI)
shading interp;
axis ij image;
colorbar
title('PHI')

subplot(4,5,11)
pcolor(xp/1000,yp/1000,pr-pf)
shading interp;
axis ij image;
colorbar
title('Ptotal-Pfluid, Pa')

subplot(4,5,12)
pcolor(xvx/1000,yvx/1000,vxf)
shading interp;
axis ij image;
colorbar
title('Vxfluid, m/s')
 

subplot(4,5,13)
pcolor(xvy/1000,yvy/1000,vyf)
shading interp;
axis ij image;
colorbar
title('VyFluid, m/s')

subplot(4,5,14)
pcolor(xvx/1000,yvx/1000,XHT)
shading interp;
axis ij image;
colorbar
title('XHT')


subplot(4,5,15)
pcolor(xp/1000,yp/1000,XHS)
shading interp;
axis ij image;
colorbar
title('XHS')

subplot(4,5,16)
pcolor(xp/1000,yp/1000, ALPHA)
shading interp;
axis ij image;
colorbar
title(' ALPHA, 1/K')


subplot(4,5,17)
pcolor(xp/1000,yp/1000, ALPHAFLUID)
shading interp;
axis ij image;
colorbar
title('ALPHA FLUID, 1/K')

subplot(4,5,18)
pcolor(xp/1000,yp/1000, ALPHASOLID)
shading interp;
axis ij image;
colorbar
title('ALPHA SOLID, 1/K')

subplot(4,5,19)
pcolor(xxm/1000,yym/1000,mtyp)
shading flat;
axis ij image;
colorbar
title('Rock type')

subplot(4,5,20)
pcolor(xxm/1000,yym/1000,mpor)
shading flat;
axis ij image;
colorbar
title('porosity')

if(1==0)

figure(5);colormap('Jet');clf
subplot(2,3,1)
pcolor(xvx/1000,yvx/1000,log10(RX))
shading interp;
axis ij image;
colorbar
title('logRX')

subplot(2,3,2)
pcolor(xvy/1000,yvy/1000,log10(RY))
shading interp;
axis ij image;
colorbar
title('logRY')

subplot(2,3,3)
pcolor(xvx/1000,yvx/1000,RHOFX)
shading interp;
axis ij image;
colorbar
title('FROFX')

subplot(2,3,4)
pcolor(xvy/1000,yvy/1000,RHOFY)
shading interp;
axis ij image;
colorbar
title('RHOFY')

subplot(2,3,5)
pcolor(xp/1000,yp/1000,log10(1./(ETAPHI+1e-30)))
shading interp;
axis ij image;
colorbar
title('logETAPHI')

subplot(2,3,6)
pcolor(xp/1000,yp/1000,log10(PHI))
shading interp;
axis ij image;
colorbar
title('log10(PHI)')

figure(6);colormap('Jet');clf
subplot(2,2,1)
pcolor(xxm,yym,mtyp)
shading flat;
axis ij image;
colorbar
title('Rock type')

subplot(2,2,2)
pcolor(xxm,yym,log10(mdis+0.001))
shading flat;
axis ij image;
colorbar
title('log10(distance), m')

subplot(2,2,3)
pcolor(xxm,yym,abs(mnum))
shading flat;
axis ij image;
colorbar
title('marker number')

subplot(2,2,4)
pcolor(xxm,yym,mpor)
shading flat;
axis ij image;
colorbar
title('porosity')

pause(0.01)
if(timestep==5 || fix(timestep/savefigstep)*savefigstep==timestep)
    namejpg    =  [nname,num2str(timestep)];
    print ('-djpeg', '-r300',namejpg);
end
end

end

if(timestep==5 || fix(timestep/savematstep)*savematstep==timestep)
    namemat    =  [nname,num2str(timestep)];
    save(namemat);
    fdata=fopen('file.txt','wt');
    fprintf(fdata,'%d',timestep);
    fclose(fdata);
end

%frame = getframe(gcf);
%writeVideo(v, frame);

end

% Close video
%close(v);