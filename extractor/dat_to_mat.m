data_name = 'BSE';

fname  = ['table_' data_name '.dat'];
outmat = ['database_' data_name '.mat'];

%% === READ TABLE ===
opts = detectImportOptions(fname,'NumHeaderLines',0);
Ttbl = readtable(fname,opts);

% New expected header order (17 columns):
%  1: T(K)
%  2: P(GPa)
%  3: H₂O_mol
%  4: Solid-H⁰(kJ/mol)
%  5: Liquid-H⁰(kJ/mol)
%  6: cptot(kJ/K/mol)          % New column for specific heat capacity
%  7: β_solid(1/GPa)
%  8: Solid-ρ(kg/m³)
%  9: Sys-ρ(kg/m³)
% 10: XH2Osolid
% 11: α_liquid(1/K)
% 12: β_liquid(1/GPa)
% 13: α_solid(1/K)
% 14: Liquid-ρ(kg/m³)
% 15: Sys-H⁰(kJ/mol)
% 16: Vtotal(cm³/mol)
% 17: X(mol_frac)

A = Ttbl{:,:};

T_K           = A(:,1);
P_GPa         = A(:,2);
H2O_mol       = A(:,3);    % hydration variable (mol or mol fraction)

Hsolid_kJmol  = A(:,4);
Hliquid_kJmol = A(:,5);

cptot_col     = A(:,6);    % New column for specific heat capacity

b_solid_col   = A(:,7);    % β_solid(1/GPa)
rhosolid_col  = A(:,8);    % Solid-ρ
rhosys_col    = A(:,9);    % Sys-ρ
XH2Osolid_col = A(:,10);   % XH2Osolid

a_liquid_col  = A(:,11);   % α_liquid(1/K)
b_liquid_col  = A(:,12);   % β_liquid(1/GPa)
a_solid_col   = A(:,13);   % α_solid(1/K)

rholiquid_col = A(:,14);   % Liquid-ρ
Hsys_kJmol    = A(:,15);   % Sys-H⁰
Vtotal_col    = A(:,16);   % Vtotal
Xmol_col      = A(:,17);   % X(mol_frac)

%% === BUILD AXES (unique sorted values) ===
tkj    = sort(uniquetol(T_K(:),    1e-9));  % Temperature axis (K)
pbi    = sort(uniquetol(P_GPa(:),  1e-9));  % Pressure axis (GPa)
wth20k = sort(uniquetol(H2O_mol(:),1e-9));  % Hydration axis (H2O_mol)

%% === COUNTS ===
tknum  = numel(tkj);    % number of temperature points
pbnum  = numel(pbi);    % number of pressure points
nwtnum = numel(wth20k); % number of hydration points

%% === STEP SIZES & MINIMA (clean + robust) ===

% Temperature step
if tknum > 1
    dT = diff(tkj);
    dT = dT(dT > eps(max(abs(tkj)))); 
    if isempty(dT)
        tkstp = 0;
    else
        dTuniq = uniquetol(dT,1e-6);
        tkstp  = min(dTuniq);
    end
else
    tkstp = 0;
end

% Pressure step
if pbnum > 1
    dP = diff(pbi);
    dP = dP(dP > eps(max(abs(pbi))));
    if isempty(dP)
        pbstp = 0;
    else
        dPuniq = uniquetol(dP,1e-6);
        pbstp  = min(dPuniq);
    end
else
    pbstp = 0;
end

% Hydration step
if nwtnum > 1
    dW = diff(wth20k);
    dW = dW(dW > eps(max(abs(wth20k))));
    if isempty(dW)
        Wth20stp = 0;
    else
        dWuniq   = uniquetol(dW,1e-6);
        Wth20stp = min(dWuniq);
    end
else
    Wth20stp = 0;
end

tkmin    = min(tkj);
pbmin    = min(pbi);
wth20min = min(wth20k);

%% === MAP EACH ROW TO (ip,it,iw) ===
[~, it] = ismember(T_K,     tkj);
[~, ip] = ismember(P_GPa,   pbi);
[~, iw] = ismember(H2O_mol, wth20k);

%% === PREALLOCATE 3D ARRAYS (np × nt × nwt) ===
rholiquid = nan(pbnum, tknum, nwtnum);
rhosolid  = nan(pbnum, tknum, nwtnum);
rhosys    = nan(pbnum, tknum, nwtnum);

Hsolid    = nan(pbnum, tknum, nwtnum);
Hliquid   = nan(pbnum, tknum, nwtnum);
Hsys      = nan(pbnum, tknum, nwtnum);

XH2Osolid = nan(pbnum, tknum, nwtnum);
Vtotal    = nan(pbnum, tknum, nwtnum);
Xmol      = nan(pbnum, tknum, nwtnum);

a_liquid  = nan(pbnum, tknum, nwtnum);
b_liquid  = nan(pbnum, tknum, nwtnum);
a_solid   = nan(pbnum, tknum, nwtnum);
b_solid   = nan(pbnum, tknum, nwtnum);

cptot     = nan(pbnum, tknum, nwtnum); % New 3D array for specific heat capacity

%% === FILL ARRAYS ===
N = size(A,1);
for k = 1:N
    if ip(k)>0 && it(k)>0 && iw(k)>0
        % densities
        rholiquid(ip(k), it(k), iw(k)) = rholiquid_col(k);
        rhosolid( ip(k), it(k), iw(k)) = rhosolid_col(k);
        rhosys(   ip(k), it(k), iw(k)) = rhosys_col(k);

        % enthalpies (kJ/mol)
        Hsolid(  ip(k), it(k), iw(k))  = Hsolid_kJmol(k);
        Hliquid( ip(k), it(k), iw(k))  = Hliquid_kJmol(k);
        Hsys(    ip(k), it(k), iw(k))  = Hsys_kJmol(k);

        % composition / volumes
        XH2Osolid(ip(k), it(k), iw(k)) = XH2Osolid_col(k);
        Vtotal(   ip(k), it(k), iw(k)) = Vtotal_col(k);
        Xmol(     ip(k), it(k), iw(k)) = Xmol_col(k);

        % thermal expansion / compressibility
        a_liquid(ip(k), it(k), iw(k))  = a_liquid_col(k);
        b_liquid(ip(k), it(k), iw(k))  = b_liquid_col(k);
        a_solid( ip(k), it(k), iw(k))  = a_solid_col(k);
        b_solid( ip(k), it(k), iw(k))  = b_solid_col(k);

        % specific heat capacity (cptot)
        cptot(ip(k), it(k), iw(k)) = cptot_col(k);
    end
end

% (Optional) Convert enthalpies to J/mol:
% Hsolid  = 1000 * Hsolid;
% Hliquid = 1000 * Hliquid;
% Hsys    = 1000 * Hsys;

%% === SAVE TO MAT FILE ===
save(outmat, ...
    'rholiquid','rhosolid','rhosys', ...
    'Hsolid','Hliquid','Hsys', ...
    'XH2Osolid','Vtotal','Xmol', ...
    'a_liquid','b_liquid','a_solid','b_solid', ...
    'cptot', ... % New field for specific heat capacity
    'tkj','pbi','wth20k', ...
    'tknum','pbnum','nwtnum', ...
    'tkstp','pbstp','Wth20stp', ...
    'tkmin','pbmin','wth20min');

disp("done");
