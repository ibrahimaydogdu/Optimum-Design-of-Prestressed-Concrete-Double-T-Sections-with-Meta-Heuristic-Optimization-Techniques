function [p_cost,penalty,fitness,sect] = Computation(sect)
%% data required
global design_Method write
global fp_c fp_ci fci fti fcs fts An_ogd nu f_pu conc_grad nu_M section Section_name...
    Section_Name Zt Zb A I yt yb h B C hf Tt Tb e tendon dp dpten n_w diam tend_pr W_p conc_pr

% Ratio of instantaneous prestress losses nu
nu = nu_M(sect(1));

% Double T beam sections 
Zt = section(sect(2), 1); Zb = section(sect(2), 2); A = section(sect(2), 3);
I = section(sect(2), 4); yt = section(sect(2), 5); yb = section(sect(2), 6);
h = section(sect(2), 7); B = section(sect(2), 8); C = section(sect(2), 9);
hf = section(sect(2), 10); Tt = section(sect(2), 11); Tb = section(sect(2), 12);

% Name of section
Section_Name = Section_name(sect(2),:);

% Tendon properties
n_w = tendon(sect(3), 1); % Number of strands
f_pu = tendon(sect(3), 2); % Specified tensile strength of prestressing steel
diam = tendon(sect(3), 3); % Nominal diameter of bar, wire, or prestressing strand 
An_ogd = tendon(sect(3), 4); % Nominal diameter of bar, wire, or prestressing strand
W_p = tendon(sect(3), 5); % Nominal weight of strands
tend_pr = tendon(sect(3), 6); % Nominal price of strands

% Distance from extreme compression fiber to centroid of longitudinal tension reinforcement 
dp = dpten(sect(4)) * h;
% Eccentricity of design load or prestressing force parallel to axis measured from centroid of section or reaction
e = dp - yt;

% Specified compressive strength of concrete
fp_c = conc_grad(sect(5), 1);

% Specified compressive strength of concrete at time of initial prestress
fp_ci = conc_grad(sect(5), 2);

% Concrete price
conc_pr = conc_grad(sect(5), 3);
% conc_pr = 220; % Ex 2

%% Stress boundaries
fci = - 0.6 * fp_ci; % Admissible compression stress at transfer (MPa)
fti = 0.25 * sqrt(fp_ci); % Admissible tensile stress at transfer(MPa)
fcs = - 0.45 * fp_c; % Admissible compression stress at the serviceability limit state (MPa)
fts = 0.5 * sqrt(fp_c); % Admissible tensile stress at the serviceability limit state (MPa)

%% Computation function's calling
if design_Method == 1     
    % Design of Beam with Variable Eccentricity Tendons
    [p_cost,penalty,fitness,sect] = VET(sect, write);
elseif design_Method == 2
    % Design of Beam with Constant Eccentricity Tendons
    [p_cost,penalty,fitness,sect] = CET(sect, write);
else
    % Design of Beam with Deflected Eccentricity Tendons
    [p_cost,penalty,fitness,sect] = DET(sect, write);
end
if p_cost<0;p_cost=1e30;fitness=1/p_cost;end
end

%% Design of Beam with Variable Eccentricity Tendons
function [p_cost,penalty,fitness,sect] = VET(sect, write)
%% Required Input variables
global L w_l w_sdl gama_p beta1 fci fti fcs fts Section_Name fp_c fp_ci gama_conc An_ogd nu f_pu...
    Zt Zb A I yt yb h B C hf Tt Tb e dp n_w diam tend_pr W_p conc_pr n_s...
    Variable_section_properties Opt_Method

%% Midspan moment for recomputed design loads
wo = (gama_conc * A * 1e-6); % (kN/m)
wd = w_sdl * (B + 2*C) * 1e-3; % (kN/m)
wl = w_l * (B + 2*C) * 1e-3; % (kN/m)

Mo = (wo * L^2) / 8; % (kN.m)
Md = (wd * L^2) / 8; % (kN.m)
Ml = (wl * L^2) / 8; % (kN.m)
Mt = 1.1 * (Mo + Md) + Ml ; % (kN.m)
Mu = 1.4 * (Mo + Md) + 1.7 * Ml ; % (kN.m)

%% Minimum required values for the section modulus
Z_top_min= ((Mo * (1 - nu) + Md + Ml) * 1e6) / (nu*fti - fcs); % (mm^3)
Z_bot_min= ((Mo * (1 - nu) + Md + Ml) * 1e6) / (fts - nu*fci); % (mm^3)

% pen_sec
if abs(Zt) >= abs(Z_top_min)
    pen_sec1 = 0;
 else
    pen_sec1 = abs(Z_top_min) / abs(Zt) - 1;
 end
if abs(Zb) >= abs(Z_bot_min)
	pen_sec2 = 0;
else
	pen_sec2 = abs(Z_bot_min) / abs(Zb) - 1;
end
pen_sec = pen_sec1 + pen_sec2;

%% Prestressing Force at Transfer Pi
% fcci
fcci = fti - (yt / h) * (fti-fci); % (MPa)

% Pi
r2 = I / A; % (mm^2)
P_i = ((nu *(1 + e * yb / r2) / (- fts + Mt * 1e6 / Zb) / A) ^-1); % (N)

% emax
e_max = (fti - fcci) * (Zt / P_i) + (Mo * 1e6 / P_i); % (mm)

% Pen_e
if abs(e) / abs(e_max) <= 1
    pen_e = 0;
else
    pen_e = e / e_max - 1;
end

%% Effective prestress force after all losses
P_e = nu * P_i; % (N)

%% Specified yield strength of prestressing steel 
f_py = 0.85 * f_pu; % (MPa)

f_py_a = 0.82 * f_py; % (MPa)

f_pu_a = 0.74 * f_pu; % (MPa)

%%  Area of prestressing steel
A_tendon = P_i / min(f_py_a, f_pu_a) ; % (mm^2)

%% Number of strands
n_s = round(A_tendon / An_ogd + 0.5);
if n_s>20
    penn = n_s/20-1;
else
    penn=0;
end

%% New tendon
An_tendon = n_s * An_ogd; % (mm^2)
Pn_i = 0.82 * 0.85 * f_pu * An_tendon; % (N)
Pn_e = nu * Pn_i; % (N)

%% Flexural control
% Stress in prestressed reinforcement at nominal strength of component
ho_p = An_tendon / ((B + 2*C) * (dp)); % (birimsiz)
f_ps = f_pu * (1 - (gama_p / beta1) * (ho_p * f_pu) / fp_c); % (MPa)

% Depth of equivalent rectangular stress block
a = An_tendon * f_ps / (0.85 * fp_c * (B + 2*C)); % (mm)

% Nominal flexural strength at section
Mn = An_tendon * f_ps * (dp - a / 2) * 1e-6; % (kN.m)
fi_Mn = 0.9 * Mn; % (kN.m)

% pen_flexural
if abs(Mu) / abs(fi_Mn) <= 1
    pen_flexural = 0;
else
    pen_flexural = abs(Mu) / abs(fi_Mn) - 1;
end

%% Cracking control
Mcr = (0.35 * sqrt(fp_c) * Zb + Pn_e * (I / A / yb + e)) * 1e-6; % (kN.m)
Mcr_k = 1.2 * Mcr; % (kN.m)

% pen_Cracking
if abs(Mcr_k) / abs(fi_Mn) <= 1
    pen_cr = 0;
else
    pen_cr = abs(Mcr_k) / abs(fi_Mn) - 1;
end

%% Shear Control
% Distance from extreme compression fiber to centroid of prestressing steel in tensionzone
d_p = max(0.8 * h, dp); % (mm)

% Web width
b_w = 12.5 * 25.4; % (mm) % cf kitabin sonunda

% Shear forces
Vo = wo * L / 2; % (kN)
D_Vu = (wd + wl) * L / 2; % (kN)

% Moments
D_Mu = (wd + wl) * L^2 / 8; % (kN.m)
D_Mcr = Mcr - Mo; % (kN.m)

% Nominal shear strength provided by concrete when diagonal cracking results from combined shear and moment
V_ci1 = 0.05 * sqrt(fp_c) * b_w * d_p + Vo * 1e3 + D_Vu / D_Mu * D_Mcr * 1e3; % (N)
V_ci2 = 1.7 * sqrt(fp_c) * b_w * d_p; % (N)
V_ci = max(V_ci1, V_ci2) * 1e-3; % (kN)

% Nominal shear strength provided by concrete when diagonal cracking results from high principal tensile stress in the web 
V_cw = (0.29 * sqrt(fp_c) + 0.3 * P_e / A) * b_w * d_p * 1e-3; % (kN)

% Nominal shear strength provided by concrete
V_c = min(V_ci, V_cw); % (kN)

% Factored shear force at section
V_u = (1.4 * (wo + wd) + 1.7 * wl) * L / 2; % (kN)

% Area of shear reinforcement within spacing s=200
A_v1 = (V_u / 0.75 - V_c) * 1e3 * 200 / 420 / d_p; % (mm^2/200)
A_v2 = A_tendon / 80 * f_pu / 420 * 200 / d_p * sqrt(d_p / b_w); % (mm^2/200) 
A_v3 = 0.062 * sqrt(fp_c) * b_w * 200 / 420; % (mm^2/200)
A_v4 = 0.35 * b_w * 200 / 420; % (mm^2/200)
A_v = max([A_v1, A_v2, A_v3, A_v4]); % (mm^2/200)

% db_t = sqrt(A_v * 4 / pi) / 2; % (mm)

% Nominal shear strength provided by shear reinforcement
V_s = A_v / 200 * 420 * d_p * 1e-3; % (kN)

% pen_Shear
fi_Vt = 0.75 * (V_c + V_s); % (kN)

if abs(V_u) / abs(fi_Vt) < 1
    pen_Shear = 0;
else
    pen_Shear = abs(V_u) / abs(fi_Vt) - 1;
end

%% ACI Allowable Stresses(Midspan)
% ***Transmission***

% Pi/Ac
fb_MA_Pi = -P_i / A; % (MPa)
ft_MA_Pi = -P_i / A; % (MPa)

% Pie/Z
fb_MA_Pie = -P_i * e / Zb; % (MPa)
ft_MA_Pie = P_i * e / Zt; % (MPa)

% Mo/z
fb_MA_Mo = Mo * 1e6 / Zb; % (MPa)
ft_MA_Mo = -Mo * 1e6 / Zt; % (MPa)

% Total
fb_MA_t = fb_MA_Pi + fb_MA_Pie + fb_MA_Mo; % (MPa)
ft_MA_t = ft_MA_Pi + ft_MA_Pie + ft_MA_Mo; % (MPa)

% Allowable Sigma
% For fb
if fb_MA_t > 0
    fb_MA_se = fti; % (MPa)
else
    fb_MA_se = fci; % (MPa)
end
% For ft
if ft_MA_t > 0
    ft_MA_se = fti; % (MPa)
else
    ft_MA_se = fci; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_MA_se) >= abs(fb_MA_t)
    pen_MA1 = 0;
else
    pen_MA1 = abs(fb_MA_t) / abs(fb_MA_se) - 1;
end
% For ft
if abs(ft_MA_se) >= abs(ft_MA_t)
	pen_MA2 = 0;
else
	pen_MA2 = abs(ft_MA_t) / abs(ft_MA_se) - 1;
end
pen_MA = pen_MA1 + pen_MA2;

% ***Effective***

% Pe/Ac
fb_ME_Pe = -P_e / A; % (MPa)
ft_ME_Pe = -P_e / A; % (MPa)

% Pee/Z
fb_ME_Pee = -P_e * e / Zb; % (MPa)
ft_ME_Pee = P_e * e / Zt; % (MPa)

% Mt/z
fb_ME_Mt = Mt * 1e6 / Zb; % (MPa)
ft_ME_Mt = -Mt * 1e6 / Zt; % (MPa)

% Total
fb_ME_t = fb_ME_Pe + fb_ME_Pee + fb_ME_Mt; % (MPa)
ft_ME_t = ft_ME_Pe + ft_ME_Pee + ft_ME_Mt; % (MPa)

% Allowable Sigma
% For fb
if fb_ME_t > 0
    fb_ME_se = fts; % (MPa)
else
    fb_ME_se = fcs; % (MPa)
end
% for ft
if ft_ME_t > 0
    ft_ME_se = fts; % (MPa)
else
    ft_ME_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_ME_se) >= abs(fb_ME_t)
    pen_ME1 = 0;
else
    pen_ME1 = abs(fb_ME_t) / abs(fb_ME_se) - 1;
end
% For ft
if abs(ft_ME_se) >= abs(ft_ME_t)
	pen_ME2 = 0;
else
	pen_ME2 = abs(ft_ME_t) / abs(ft_ME_se) - 1;
end
pen_ME = pen_ME1 + pen_ME2;

%% Aci allowable stress(Support)
% ***Transmission***

% Total
fb_SA_t = -Pn_i / A; % (MPa)
ft_SA_t = fb_SA_t; % (MPa)

% Allowable Sigma
% for fb
if fb_SA_t > 0
    fb_SA_se = 2 * fti; % (MPa)
else
    fb_SA_se = fci; % (MPa)
end
% For ft
if ft_SA_t > 0
    ft_SA_se = 2 * fti; % (MPa)
else
    ft_SA_se = fci; % (MPa)
end

% pen_ACI_S
% For fb
if abs(fb_SA_se) >= abs(fb_SA_t)
    pen_SA1 = 0;
else
    pen_SA1 = abs(fb_SA_t) / abs(fb_SA_se) - 1;
end
% For ft
if abs(ft_SA_se) >= abs(ft_SA_t)
	pen_SA2 = 0;
else
	pen_SA2 = abs(ft_SA_t) / abs(ft_SA_se) - 1;
end
pen_SA = pen_SA1 + pen_SA2;

% ***Effective***

% Total
fb_SE_t = -Pn_e / A; % (MPa)
ft_SE_t = fb_SE_t; % (MPa)

% Allowable Sigma
% For fb
if fb_SE_t > 0
    fb_SE_se = fts; % (MPa)
else
    fb_SE_se = fcs; % (MPa)
end
% For ft
if ft_SE_t > 0
    ft_SE_se = fts; % (MPa)
else
    ft_SE_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_SE_se) >= abs(fb_SE_t)
    pen_SE1 = 0;
else
    pen_SE1 = abs(fb_SE_t) / abs(fb_SE_se) - 1;
end
% For ft
if abs(ft_SE_se) >= abs(ft_SE_t)
	pen_SE2 = 0;
else
	pen_SE2 = abs(ft_SE_t) / abs(ft_SE_se) - 1;
end
pen_SE = pen_SE1 + pen_SE2;

%% Costs 
% Cost of concrete
conc_cost = conc_pr * (A * 1e-6);

% Cost of tendon
tendon_cost = n_s * n_w * tend_pr * W_p;

% Cost of Reinforcement Steel
   % For a price of reinforcement steel 250 Euro/m^3, we have:
reinf_cost = 250 * A_v * 1e-6;

% Total Cost
Cost = conc_cost + tendon_cost + reinf_cost ; % Cost per metre

%% X Evaluate the Design
pen_str = pen_MA + pen_ME + pen_SA + pen_SE;
pen_UM = pen_flexural + pen_cr;
penalty = pen_sec + pen_e + pen_UM + pen_Shear + pen_str + penn;
p_cost = Cost * (1 + penalty)^2;
fitness = 1 / p_cost;

%% Displaying
if write == 1
    disp('---------------------------------------------------------')
    if Opt_Method == 1
        disp('<strong>********** ABC Design of Beam with Variable Eccentricity Tendons **********</strong>')
    elseif Opt_Method == 2
        disp('<strong>********** BSO Design of Beam with Variable Eccentricity Tendons **********</strong>')
    else
        disp('<strong>********** SSO Design of Beam with Variable Eccentricity Tendons **********</strong>')
    end
    disp('<strong>Section Properties</strong>')
    fprintf('Selected Section Name : <strong>%s\n</strong>', Section_Name(:))
    fprintf('Zt = %7.4f mm^3\n', Zt)
    fprintf('Zb = %7.4f mm^3\n', Zb)
    fprintf('A  = %7.4f mm^2\n', A)
    fprintf('I  = %7.4f mm^4\n', I)
    fprintf('yt = %7.4f mm\n', yt)
    fprintf('yb = %7.4f mm\n', yb)
    fprintf('h  = %7.4f mm\n', h)
    fprintf('B  = %7.4f mm\n', B)
    fprintf('C  = %7.4f mm\n', C)
    fprintf('hf = %7.4f mm\n', hf)
    fprintf('Tt = %7.4f mm\n', Tt)
    fprintf('Tb = %7.4f mm\n', Tb)
    fprintf('Position of prestress tendon = %7.4f mm\n', dp)
    fprintf('                Eccentricity = %7.4f mm\n', e)
    fprintf('  Initial Prestressing Force = %7.4f kN\n', P_i*1e-3)
    fprintf('  Effective prestress force  = %7.4f kN\n', P_e*1e-3)
    disp('---------------------------------------------------------')
    disp('<strong>Tendon Properties</strong>')
    fprintf('  nwire = %7.4f \n', n_w)
    fprintf('    fpu = %7.4f MPa\n', f_pu)
    fprintf('     db = %7.4f mm^2\n', diam)
    fprintf('ntendon = %7.4f \n', n_s)
    fprintf('     Ap = %7.4f mm^2\n', An_ogd)
    fprintf('Atendon = %7.4f mm^2\n', An_tendon)
    fprintf('     Wp = %7.4f kg/m\n', W_p)

    disp('---------------------------------------------------------')
    disp('  <strong>I. Minimum required values for the section modulus</strong>')
    if Zt >= abs(Z_top_min) && Zb >= abs(Z_bot_min)
        fprintf( 'Zt = %3.4g mm > Z_top_min = %3.4g mm\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm > Z_bot_min = %3.4g mm\n', Zb, abs(Z_bot_min))
        disp('<strong>OK</strong>');
    else    
        fprintf( 'Zt = %3.4g mm >= Z_top_min = %3.4g mm olmasi lazim\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm >= Z_bot_min = %3.4g mm olmasi lazim\n', Zb, abs(Z_bot_min))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>II. Flexural Control</strong>')
    if Mu / fi_Mn <= 1
        fprintf('Mu = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mu = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>NOT OK</strong>')
    end

    disp(' ')
    disp('  <strong>III. Cracking Control</strong>')
    if Mcr_k / fi_Mn <= 1
        fprintf('Mcr_k = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mcr_k = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>IV. Shear Control</strong>')
    if  V_u / fi_Vt <= 1
        fprintf('V_u = %3.4g kN <= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>OK</strong>');
    else
        fprintf('V_u = %3.4g kN >= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>V. Aci Allowable Stresses</strong>')
    disp(' <strong>V.1 Midspan</strong>')
    disp('<strong>V.1.1 Transmission</strong>')

    % For fb
    if  abs(fb_MA_t) / abs(fb_MA_se) <= 1
        fprintf('abs(fb_MA_t) = %3.4g MPa <= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_MA_t) = %3.4g MPa >= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if  abs(ft_MA_t) / abs(ft_MA_se) <= 1
        fprintf('abs(ft_MA_t) = %3.4g MPa <= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_MA_t) = %3.4g MPa >= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('<strong>V.1.2 Effective</strong>')
    % For fb
    if  abs(fb_ME_t) / abs(fb_ME_se) <= 1
        fprintf('abs(fb_ME_t) = %3.4g MPa <= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_ME_t) = %3.4g MPa >= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_ME_t) / abs(ft_ME_se) <= 1
        fprintf('abs(ft_ME_t) = %3.4g MPa <= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_ME_t) = %3.4g MPa >= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp(' <strong>V.2 Support</strong>')
    disp('<strong>V.2.1 Transmission</strong>')
    % For fb
    if  abs(fb_SA_t) / abs(fb_SA_se) <= 1
        fprintf('abs(fb_SA_t) = %3.4g MPa <= abs(fb_SA_se) = %3.4g MPa\n', abs(fb_SA_t), abs(fb_SA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_SA_t) = %3.4g MPa >= abs(fb_SA_se) = %3.4g MPa\n', abs(fb_SA_t), abs(fb_SA_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if  abs(ft_SA_t) / abs(ft_SA_se) <= 1
        fprintf('abs(ft_SA_t) = %3.4g MPa <= abs(ft_SA_se) = %3.4g MPa\n', abs(ft_SA_t), abs(ft_SA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_SA_t) = %3.4g MPa >= abs(ft_SA_se) = %3.4g MPa\n', abs(ft_SA_t), abs(ft_SA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('<strong>V.2.2 Effective</strong>')
    % For fb
    if  abs(fb_SE_t) / abs(fb_SE_se) <= 1
        fprintf('abs(fb_SE_t) = %3.4g MPa <= abs(fb_SE_se) = %3.4g MPa\n', abs(fb_SE_t), abs(fb_SE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_SE_t) = %3.4g MPa >= abs(fb_SE_se) = %3.4g MPa\n', abs(fb_SE_t), abs(fb_SE_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_SE_t) / abs(ft_SE_se) <= 1
        fprintf('abs(ft_SE_t) = %3.4g MPa <= abs(ft_SE_se) = %3.4g MPa\n', abs(ft_SE_t), abs(ft_SE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_SE_t) = %3.4g MPa >= abs(ft_SE_se) = %3.4g MPa\n', abs(ft_SE_t), abs(ft_SE_se))
        disp('<strong>NOT OK</strong>');
    end

    disp('---------------------------------------------------------');
    disp('<strong>Penalty</strong>')
    fprintf('      Penalty from Section = %7.4f \n',pen_sec)
    fprintf(' Penalty from Eccentricity = %7.4f \n', pen_e)
    fprintf(' Penalty from Moment_Check = %7.4f \n', pen_UM)
    fprintf('  Penalty from Shear_Check = %7.4f \n', pen_Shear)
    fprintf('Penalty from Tendon stress = %7.4f \n', pen_str)

    %% Table
    section_properties_names = {'fcp', 'fcip', 'Gama conc', 'Ratio nu', 'Zt', 'Zb', 'A', 'I', 'yt', 'yb', 'h', 'B', 'C', 'hf', 'Tt', 'Tb','Zt_min', 'Zb_min'}';
    section_properties_values = [fp_c, fp_ci, gama_conc, nu, Zt, Zb, A, I, yt, yb, h, B, C, hf, Tt, Tb Z_top_min Z_bot_min]';
    section_properties_units = {'Mpa', 'Mpa', 'kN/m^3', '-', 'mm^3', 'mm^3', 'mm^2', 'mm^4', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm^3', 'mm^3'}';
    
    Tendon_properties_names = {'dp', 'e', 'Pi', 'Pe', 'nwire', 'fpu', 'db', 'n_tendon', 'Ap', 'A_tendon', 'Wp'}';
    Tendon_properties_values = [dp e P_i*1e-3 P_e*1e-3 n_w f_pu diam n_s An_ogd An_tendon W_p]';
    Tendon_properties_units = {'mm', 'mm', 'kN', 'kN', '-', 'Mpa', 'mm^2', '-', 'mm^2', 'mm^2', 'kg/m'}';
    
    Moment_kontrol_names = {'Mu', 'fi_Mn', 'Mcr', '1.2Mcr', 'Vu', 'fi_Vt'}';
    Moment_kontrol_values = [Mu fi_Mn Mcr Mcr_k V_u fi_Vt]';
    Moment_kontrol_units = {'kN.m', 'kN.m', 'kN.m', 'kN.m', 'kN', 'kN'}';
    
    Cost_name = {'Cost'};
    Cost_value = Cost;
    Cost_unit = {'$'};
    Rownames = [section_properties_names; Tendon_properties_names; Moment_kontrol_names; Cost_name];
    Variable = [section_properties_values; Tendon_properties_values; Moment_kontrol_values; Cost_value];
    Units = categorical([section_properties_units; Tendon_properties_units; Moment_kontrol_units; Cost_unit]);
    Variable_section_properties = table(Variable, Units, 'RowNames',Rownames);
       
    %% Limit Zone For Tendon Centroid
    %**Midspan**%

    % e
    M_e = (-200:5:900)';
    M_size = size(M_e,1);

    % 1/Pi>fti
    M_1_bo_Pi_bu_fti(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fti(i,1) = (-1 + M_e(i) * yt / I * A) / (fti + Mo * 1e6 / Zt) / A; 
    end

    % 1/Pi>fci
    M_1_bo_Pi_bu_fci(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fci(i,1) = (1 + M_e(i) * yb / I * A) / (-fci + Mo * 1e6 / Zb) / A; 
    end

    % 1/Pi<fts
    M_1_bo_Pi_ku_fts(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + M_e(i) * yb / I * A) / (-fts + Mt * 1e6 / Zb) / A; 
    end

    % 1/Pi<fcs
    M_1_bo_Pi_ku_fcs(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + M_e(i) * yt / I * A) / (fcs + Mt * 1e6 / Zt) / A; 
    end

    % Pi
    M_Pi(M_size,4)=0;
    for i = 1:M_size
        M_Pi(i,1) = 1 / M_1_bo_Pi_bu_fti(i);
        M_Pi(i,2) = 1 / M_1_bo_Pi_bu_fci(i);
        M_Pi(i,3) = 1 / M_1_bo_Pi_ku_fts(i);
        M_Pi(i,4) = 1 / M_1_bo_Pi_ku_fcs(i);
    end

    M_Max = max([M_1_bo_Pi_bu_fti; M_1_bo_Pi_bu_fci; M_1_bo_Pi_ku_fts; M_1_bo_Pi_ku_fcs]);

    %***Graphic 1***%
    figure(1)
        plot(M_1_bo_Pi_bu_fti, M_e)
        title('Design of Beam with Variable Eccentricity Tendons - Midspan')
        ylabel('Eccentricity, mm')
        ylim([M_e(1) - 200 M_e(end) + 200])
        axis ij
        xlim([0 M_Max])
        grid on
        hold on
        plot(M_1_bo_Pi_bu_fci, M_e)
        plot(M_1_bo_Pi_ku_fts, M_e)
        plot(M_1_bo_Pi_ku_fcs, M_e)
        plot([0 M_Max], [e e])
        plot(1/P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off              
        
    %**Support**%
    % e
    S_e = (-200:5:900)';
    S_size = size(S_e,1);

    % 1/Pi>fti
    S_1_bo_Pi_bu_fti(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_bu_fti(i,1) = (-1 + M_e(i) * yt / I * A) / (2 * fti * A);
    end

    % 1/Pi>fci
    S_1_bo_Pi_bu_fci(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_bu_fci(i,1) = (1 + M_e(i) * yb / I * A) / (-fci * A);
    end

    % 1/Pi<fts
    S_1_bo_Pi_ku_fts(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + M_e(i) * yb / I * A) / (-fts * A);
    end

    % 1/Pi<fts
    S_1_bo_Pi_ku_fcs(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + M_e(i) * yt / I * A) / (fcs * A);
    end

    % Pi
    S_Pi(S_size,4)=0;
    for i = 1:S_size
        S_Pi(i,1) = 1 / S_1_bo_Pi_bu_fti(i);
        S_Pi(i,2) = 1 / S_1_bo_Pi_bu_fci(i);
        S_Pi(i,3) = 1 / S_1_bo_Pi_ku_fts(i);
        S_Pi(i,4) = 1 / S_1_bo_Pi_ku_fcs(i);
    end

    S_Max = max([S_1_bo_Pi_bu_fti; S_1_bo_Pi_bu_fci; S_1_bo_Pi_ku_fts; S_1_bo_Pi_ku_fcs]);

    %***Grafic 2***%
    figure(2)
        plot(S_1_bo_Pi_bu_fti, S_e)
        title('Design of Beam with Variable Eccentricity Tendons - Support')
        ylabel('Eccentricity, mm')
        ylim([S_e(1) - 200 S_e(end) + 200])
        axis ij
        xlim([0 S_Max])
        grid on
        hold on
        plot(S_1_bo_Pi_bu_fci, S_e)
        plot(S_1_bo_Pi_ku_fts, S_e)
        plot(S_1_bo_Pi_ku_fcs, S_e)
        plot([0 S_Max], [e e])
        plot(1/P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off              

    %% Saving
    % Table
    if Opt_Method == 1
        save('ABC_Variable.mat','Variable_section_properties')
    elseif Opt_Method == 2
        save('BSO_Variable.mat','Variable_section_properties')
    else
        save('SSO_Variable.mat','Variable_section_properties')
    end

    % figure 1
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Variable Eccentricity Tendons - Midspan.png')
            delete('ABC Design of Beam with Variable Eccentricity Tendons - Midspan.png')
        end
        saveas(figure(1),'ABC Design of Beam with Variable Eccentricity Tendons - Midspan.png')
        elseif Opt_Method == 2
            if isfile('BSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
                delete('BSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'BSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
        else
            if isfile('SSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
                delete('SSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'SSO Design of Beam with Variable Eccentricity Tendons - Midspan.png')
    end

    % figure 2
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Variable Eccentricity Tendons - Support.png')
            delete('ABC Design of Beam with Variable Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'ABC Design of Beam with Variable Eccentricity Tendons - Support.png')
    elseif Opt_Method == 2
        if isfile('BSO Design of Beam with Variable Eccentricity Tendons - Support.png')
            delete('BSO Design of Beam with Variable Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'BSO Design of Beam with Variable Eccentricity Tendons - Support.png')
    else
        if isfile('SSO Design of Beam with Variable Eccentricity Tendons - Support.png')
            delete('SSO Design of Beam with Variable Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'SSO Design of Beam with Variable Eccentricity Tendons - Support.png')
    end
end

end

%% Design of Beam with Constant Eccentricity Tendons
function [p_cost,penalty,fitness,sect] = CET(sect, write)
%% Required Input variables
global L w_l w_sdl gama_p beta1 fci fti fcs fts Section_Name fp_c fp_ci gama_conc An_ogd nu f_pu...
    Zt Zb A I yt yb h B C hf Tt Tb e dp n_w diam tend_pr W_p conc_pr n_s...
    Constant_section_properties Opt_Method

%% Midspan moment for recomputed design loads
wo = (gama_conc * A * 1e-6); % (kN/m)
wd = w_sdl * (B + 2*C) * 1e-3; % (kN/m)
wl = w_l * (B + 2*C) * 1e-3; % (kN/m)

Mo = (wo * L^2) / 8; % (kN.m)
Md = (wd * L^2) / 8; % (kN.m)
Ml = (wl * L^2) / 8; % (kN.m)
Mt = Mo + Md + Ml; % (kN.m) %~ DDM
Mu = 1.4 * (Mo + Md) + 1.7 * Ml ; % (kN.m)

%% Minimum required values for the section modulus
Z_top_min= Mt * 1e6 / (nu*fti - fcs); % (mm^3) %~ DDM
Z_bot_min= Mt * 1e6 / (fts - nu*fci); % (mm^3) %~ DDM

if abs(Zt)>=abs(Z_top_min)
    pen_sec1 = 0;
 else
    pen_sec1 = abs(Z_top_min) / abs(Zt) - 1;
 end
if abs(Zb)>=abs(Z_bot_min)
	pen_sec2=0;
else
	pen_sec2 = abs(Z_bot_min) / abs(Zb) - 1;
end
pen_sec = pen_sec1 + pen_sec2;

%% Prestressing Force at Transfer Pi
% fcci
fcci = fti - (yt / h) * (fti-fci); % (MPa)

% Pi
r2 = I / A; % (mm^2)
P_i = ((nu *(1 + e * yb / r2) / (- fts + Mt * 1e6 / Zb) / A) ^-1); % (N)

% emax
e_max = (fti - fcci) * (Zt / P_i); % (mm)

% Pen_e
if abs(e) / abs(e_max) <= 1
    pen_e = 0;
else
    pen_e = abs(e) / abs(e_max) - 1;
end

%% Effective prestress force after all losses
P_e = nu * P_i; % (N)

%% Specified yield strength of prestressing steel
f_py = 0.85 * f_pu; % (MPa)

f_py_a = 0.82 * f_py; % (MPa)

f_pu_a = 0.74 * f_pu; % (MPa)

%% Area of prestressing steel
A_tendon = P_i / min(f_py_a, f_pu_a); % (mm^2)

%% Number of strands
n_s = round(A_tendon / An_ogd + 0.5);
if n_s>20
    penn = n_s/20-1;
else
    penn=0;
end

%% New tendon
An_tendon = n_s * An_ogd; % (mm^2)
Pn_i = 0.82 * 0.85 * f_pu * An_tendon; % (N)
Pn_e = Pn_i * nu; % (N)

%% Flexural control
% Stress in prestressed reinforcement at nominal strength of component
ho_p = An_tendon / ((B + 2*C) * (dp)); % (birimsiz)
f_ps = f_pu * (1 - (gama_p / beta1) * (ho_p * f_pu) / fp_c); % (MPa)

% Depth of equivalent rectangular stress block
a = An_tendon * f_ps / (0.85 * fp_c * (B + 2*C)); % (mm)

% Nominal flexural strength at section
Mn = An_tendon * f_ps * (dp - a / 2) * 1e-6; % (kN.m)
fi_Mn = 0.9 * Mn; % (kN.m)

% pen_flexural
if abs(Mu) / abs(fi_Mn) <= 1
    pen_flexural = 0;
else
    pen_flexural = abs(Mu) / abs(fi_Mn) - 1;
end

%% Cracking control
Mcr = (0.35 * sqrt(fp_c) * Zb + Pn_e * (I / A / yb + e)) * 1e-6; % (kN.m)
Mcr_k = 1.2 * Mcr; % (kN.m)

% pen_Cracking
if abs(Mcr_k) / abs(fi_Mn) <= 1
    pen_cr = 0;
else
    pen_cr = abs(Mcr_k) / abs(fi_Mn) - 1;
end

%% Shear Control
% Distance from extreme compression fiber to centroid of prestressing steel in tensionzone
d_p = max(0.8 * h, dp); % (mm)

% Web width
b_w = 12.5 * 25.4; % (mm) % cf kitabin sonunda

% Shear forces
Vo = wo * L / 2; % (kN)
D_Vu = (wd + wl) * L / 2; % (kN)

% Moments
D_Mu = (wd + wl) * L^2 / 8; % (kN.m)
D_Mcr = Mcr - Mo; % (kN.m)

% Nominal shear strength provided by concrete when diagonal cracking results from combined shear and moment
V_ci1 = 0.05 * sqrt(fp_c) * b_w * d_p + Vo * 1e3 + D_Vu / D_Mu * D_Mcr * 1e3; % (N)
V_ci2 = 1.7 * sqrt(fp_c) * b_w * d_p; % (N)
V_ci = max(V_ci1, V_ci2) * 1e-3; % (kN)

% Nominal shear strength provided by concrete when diagonal cracking results from high principal tensile stress in the web 
V_cw = (0.29 * sqrt(fp_c) + 0.3 * P_e / A) * b_w * d_p * 1e-3; % (kN)

% Nominal shear strength provided by concrete
V_c = min(V_ci, V_cw); % (kN)

% Factored shear force at section
V_u = (1.4 * (wo + wd) + 1.7 * wl) * L / 2; % (kN)

% Area of shear reinforcement within spacing s=200
A_v1 = (V_u / 0.75 - V_c) * 1e3 * 200 / 420 / d_p; % (mm^2/200)
A_v2 = A_tendon / 80 * f_pu / 420 * 200 / d_p * sqrt(d_p / b_w); % (mm^2/200) 
A_v3 = 0.062 * sqrt(fp_c) * b_w * 200 / 420; % (mm^2/200)
A_v4 = 0.35 * b_w * 200 / 420; % (mm^2/200)
A_v = max([A_v1, A_v2, A_v3, A_v4]); % (mm^2/200)

% db_t = sqrt(A_v * 4 / pi) / 2; % (mm)

% Nominal shear strength provided by shear reinforcement
V_s = A_v / 200 * 420 * d_p * 1e-3; % (kN)

% pen_Shear
fi_Vt = 0.75 * (V_c + V_s); % (kN)

if abs(V_u) / abs(fi_Vt) <= 1
    pen_Shear = 0;
else
    pen_Shear = abs(V_u)/abs(fi_Vt)-1;
end

%% ACI Allowable Stresses(Midspan)
% ***Transmission***

% Pi/Ac
fb_MA_Pi = -P_i / A; % (MPa)
ft_MA_Pi = -P_i / A; % (MPa)

% Pie/Z
fb_MA_Pie = -P_i * e / Zb; % (MPa)
ft_MA_Pie = P_i * e / Zt; % (MPa)

% Mo/z
fb_MA_Mo = Mo * 1e6 / Zb; % (MPa)
ft_MA_Mo = -Mo * 1e6 / Zt; % (MPa)

% Total
fb_MA_t = fb_MA_Pi + fb_MA_Pie + fb_MA_Mo; % (MPa)
ft_MA_t = ft_MA_Pi + ft_MA_Pie + ft_MA_Mo; % (MPa)

% Allowable Sigma
% For fb
if fb_MA_t > 0
    fb_MA_se = fti; % (MPa)
else
    fb_MA_se = fci; % (MPa)
end
% For ft
if ft_MA_t > 0
    ft_MA_se = fti; % (MPa)
else
    ft_MA_se = fci; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_MA_se) >= abs(fb_MA_t)
    pen_MA1 = 0;
else
    pen_MA1 = abs(fb_MA_t) / abs(fb_MA_se) - 1;
end
% For ft
if abs(ft_MA_se) >= abs(ft_MA_t)
	pen_MA2 = 0;
else
	pen_MA2 = abs(ft_MA_t) / abs(ft_MA_se) - 1;
end
pen_MA = pen_MA1 + pen_MA2;

    
% ***Effective***

% Pe/Ac
fb_ME_Pe = -P_e / A; % (MPa)
ft_ME_Pe = -P_e / A; % (MPa)

% Pee/Z
fb_ME_Pee = -P_e * e / Zb; % (MPa)
ft_ME_Pee = P_e * e / Zt; % (MPa)

% Mt/z
fb_ME_Mt = Mt * 1e6 / Zb; % (MPa)
ft_ME_Mt = -Mt * 1e6 / Zt; % (MPa)

% Total
fb_ME_t = fb_ME_Pe + fb_ME_Pee + fb_ME_Mt; % (MPa)
ft_ME_t = ft_ME_Pe + ft_ME_Pee + ft_ME_Mt; % (MPa)

% Allowable Sigma
% For fb
if fb_ME_t > 0
    fb_ME_se = fts; % (MPa)
else
    fb_ME_se = fcs; % (MPa)
end
% for ft
if ft_ME_t > 0
    ft_ME_se = fts; % (MPa)
else
    ft_ME_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_ME_se) >= abs(fb_ME_t)
    pen_ME1=0;
else
    pen_ME1 = abs(fb_ME_t) / abs(fb_ME_se) - 1;
end
% For ft
if abs(ft_ME_se) >= abs(ft_ME_t)
	pen_ME2 = 0;
else
	pen_ME2 = abs(ft_ME_t) / abs(ft_ME_se) - 1;
end
pen_ME = pen_ME1 + pen_ME2;

%% Aci allowable stress(Support) %~ DDM
% ***Transmission***

% Pi/Ac
fb_SA_Pi = -P_i / A; % (MPa)
ft_SA_Pi = -P_i / A; % (MPa)

% Pie/Z
fb_SA_Pie = -P_i * e / Zb; % (MPa)
ft_SA_Pie = P_i * e / Zt; % (MPa)

% Total
fb_SA_t = fb_SA_Pi + fb_SA_Pie; % (MPa)
ft_SA_t = ft_SA_Pi + ft_SA_Pie; % (MPa)

% Allowable Sigma
% for fb
if fb_SA_t > 0
    fb_SA_se = 2 * fti; % (MPa)
else
    fb_SA_se = fci; % (MPa)
end
% For ft
if ft_SA_t > 0
    ft_SA_se = 2 * fti; % (MPa)
else
    ft_SA_se = fci; % (MPa)
end

% pen_ACI_S
% For fb
if abs(fb_SA_se) >= abs(fb_SA_t)
    pen_SA1=0;
else
    pen_SA1 = abs(fb_SA_t) / abs(fb_SA_se) - 1;
end
% For ft
if abs(ft_SA_se) >= abs(ft_SA_t)
	pen_SA2 = 0;
else
	pen_SA2 = abs(ft_SA_t) / abs(ft_SA_se) - 1;
end
pen_SA = pen_SA1 + pen_SA2;

% ***Effective***

% Pe/Ac
fb_SE_Pe = -P_e / A; % (MPa)
ft_SE_Pe = -P_e / A; % (MPa)

% Pee/Z
fb_SE_Pee = -P_e * e / Zb; % (MPa)
ft_SE_Pee = P_e * e / Zt; % (MPa)

% Total
fb_SE_t = fb_SE_Pe + fb_SE_Pee; % (MPa)
ft_SE_t = ft_SE_Pe + ft_SE_Pee; % (MPa)

% Allowable Sigma
% For fb
if fb_SE_t > 0
    fb_SE_se = fts; % (MPa)
else
    fb_SE_se = fcs; % (MPa)
end
% For ft
if ft_SE_t > 0
    ft_SE_se = fts; % (MPa)
else
    ft_SE_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_SE_se) >= abs(fb_SE_t)
    pen_SE1 = 0;
else
    pen_SE1 = abs(fb_SE_t) / abs(fb_SE_se) - 1;
end
% For ft
if abs(ft_SE_se) >= abs(ft_SE_t)
	pen_SE2 = 0;
else
	pen_SE2 = abs(ft_SE_t) / abs(ft_SE_se) - 1;
end
pen_SE = pen_SE1 + pen_SE2;

%% Costs 
% Cost of concrete
conc_cost = conc_pr * (A * 1e-6);

% Cost of tendon
tendon_cost = n_s * n_w * tend_pr * W_p;

% Cost of Reinforcement Steel
   % For a price of reinforcement steel 250 Euro/m^3, we have:
reinf_cost = 250 * A_v * 1e-6;

% Total Cost
Cost = conc_cost + tendon_cost + reinf_cost ; % Cost per metre

%% X Evaluate the Design

pen_str = pen_MA + pen_ME + pen_SA + pen_SE;
pen_UM = pen_flexural + pen_cr;
penalty = pen_sec + pen_e + pen_UM + pen_Shear + pen_str + penn;
p_cost = Cost * (1 + penalty)^2;
fitness = 1 / p_cost;

%% Displaying
if write == 1
    disp('---------------------------------------------------------')
    if Opt_Method == 1
        disp('<strong>********** ABC Design of Beam with Constant Eccentricity Tendons **********</strong>')
    elseif Opt_Method == 2
        disp('<strong>********** BSO Design of Beam with Constant Eccentricity Tendons **********</strong>')
    else
        disp('<strong>********** SSO Design of Beam with Constant Eccentricity Tendons **********</strong>')
    end
    disp('<strong>Section Properties</strong>')
    fprintf('Selected Section Name : <strong>%s\n</strong>', Section_Name(:))
    fprintf('Zt = %7.4f mm^3\n', Zt)
    fprintf('Zb = %7.4f mm^3\n', Zb)
    fprintf('A  = %7.4f mm^2\n', A)
    fprintf('I  = %7.4f mm^4\n', I)
    fprintf('yt = %7.4f mm\n', yt)
    fprintf('yb = %7.4f mm\n', yb)
    fprintf('h  = %7.4f mm\n', h)
    fprintf('B  = %7.4f mm\n', B)
    fprintf('C  = %7.4f mm\n', C)
    fprintf('hf = %7.4f mm\n', hf)
    fprintf('Tt = %7.4f mm\n', Tt)
    fprintf('Tb = %7.4f mm\n', Tb)
    fprintf('Position of prestress tendon = %7.4f mm\n', dp)
    fprintf('                Eccentricity = %7.4f mm\n', e)
    fprintf('  Initial Prestressing Force = %7.4f kN\n', P_i*1e-3)
    fprintf('  Effective prestress force  = %7.4f kN\n', P_e*1e-3)
    
    disp('---------------------------------------------------------')
    disp('<strong>Tendon Properties</strong>')
    fprintf('  nwire = %7.4f \n', n_w)
    fprintf('    fpu = %7.4f MPa\n', f_pu)
    fprintf('     db = %7.4f mm^2\n', diam)
    fprintf('ntendon = %7.4f \n', n_s)
    fprintf('     Ap = %7.4f mm^2\n', An_ogd)
    fprintf('Atendon = %7.4f mm^2\n', An_tendon)
    fprintf('     Wp = %7.4f kg/m\n', W_p)


    disp('---------------------------------------------------------')
    disp('  <strong>I. Minimum required values for the section modulus</strong>')
    if Zt >= abs(Z_top_min) && Zb >= abs(Z_bot_min)
        fprintf( 'Zt = %3.4g mm > Z_top_min = %3.4g mm\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm > Z_bot_min = %3.4g mm\n', Zb, abs(Z_bot_min))
        disp('<strong>OK</strong>');
    else    
        fprintf( 'Zt = %3.4g mm >= Z_top_min = %3.4g mm olmasi lazim\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm >= Z_bot_min = %3.4g mm olmasi lazim\n', Zb, abs(Z_bot_min))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>II. Flexural Control</strong>')
    if Mu / fi_Mn <= 1
        fprintf('Mu = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mu = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>NOT OK</strong>')
    end

    disp(' ')
    disp('  <strong>III. Cracking Control</strong>')
    if Mcr_k / fi_Mn <= 1
        fprintf('Mcr_k = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mcr_k = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>NOT OK</strong>');
    end
    
    disp(' ')
    disp('  <strong>IV. Shear Control</strong>')
    if  V_u / fi_Vt <= 1
        fprintf('V_u = %3.4g kN <= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>OK</strong>');
    else
        fprintf('V_u = %3.4g kN >= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>V. Aci Allowable Stresses</strong>')
    disp(' <strong>V.1 Midspan</strong>')
    disp('<strong>V.1.1 Transmission</strong>')

    % For fb
    if  abs(fb_MA_t) / abs(fb_MA_se) <= 1
        fprintf('abs(fb_MA_t) = %3.4g MPa <= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_MA_t) = %3.4g MPa >= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if  abs(ft_MA_t) / abs(ft_MA_se) <= 1
        fprintf('abs(ft_MA_t) = %3.4g MPa <= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_MA_t) = %3.4g MPa >= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('<strong>V.1.2 Effective</strong>')
    % For fb
    if  abs(fb_ME_t) / abs(fb_ME_se) <= 1
        fprintf('abs(fb_ME_t) = %3.4g MPa <= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_ME_t) = %3.4g MPa >= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_ME_t) / abs(ft_ME_se) <= 1
        fprintf('abs(ft_ME_t) = %3.4g MPa <= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_ME_t) = %3.4g MPa >= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp(' <strong>V.2 Support</strong>')
    disp('<strong>V.2.1 Transmission</strong>')
    % For fb
    if  abs(fb_SA_t) / abs(fb_SA_se) <= 1
        fprintf('abs(fb_SA_t) = %3.4g MPa <= abs(fb_SA_se) = %3.4g MPa\n', abs(fb_SA_t), abs(fb_SA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_SA_t) = %3.4g MPa >= abs(fb_SA_se) = %3.4g MPa\n', abs(fb_SA_t), abs(fb_SA_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if  abs(ft_SA_t) / abs(ft_SA_se) <= 1
        fprintf('abs(ft_SA_t) = %3.4g MPa <= abs(ft_SA_se) = %3.4g MPa\n', abs(ft_SA_t), abs(ft_SA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_SA_t) = %3.4g MPa >= abs(ft_SA_se) = %3.4g MPa\n', abs(ft_SA_t), abs(ft_SA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('<strong>V.2.2 Effective</strong>')
    % For fb
    if  abs(fb_SE_t) / abs(fb_SE_se) <= 1
        fprintf('abs(fb_SE_t) = %3.4g MPa <= abs(fb_SE_se) = %3.4g MPa\n', abs(fb_SE_t), abs(fb_SE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_SE_t) = %3.4g MPa >= abs(fb_SE_se) = %3.4g MPa\n', abs(fb_SE_t), abs(fb_SE_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_SE_t) / abs(ft_SE_se) <= 1
        fprintf('abs(ft_SE_t) = %3.4g MPa <= abs(ft_SE_se) = %3.4g MPa\n', abs(ft_SE_t), abs(ft_SE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_SE_t) = %3.4g MPa >= abs(ft_SE_se) = %3.4g MPa\n', abs(ft_SE_t), abs(ft_SE_se))
        disp('<strong>NOT OK</strong>');
    end

    disp('---------------------------------------------------------');
    disp('<strong>Penalty</strong>')
    fprintf('      Penalty from Section = %7.4f \n',pen_sec)
    fprintf(' Penalty from Eccentricity = %7.4f \n', pen_e)
    fprintf(' Penalty from Moment_Check = %7.4f \n', pen_UM)
    fprintf('  Penalty from Shear_Check = %7.4f \n', pen_Shear)
    fprintf('Penalty from Tendon stress = %7.4f \n', pen_str)

    %% Table
    section_properties_names = {'fcp', 'fcip', 'Gama conc', 'Ratio nu', 'Zt', 'Zb', 'A', 'I', 'yt', 'yb', 'h', 'B', 'C', 'hf', 'Tt', 'Tb','Zt_min', 'Zb_min'}';
    section_properties_values = [fp_c, fp_ci, gama_conc, nu, Zt, Zb, A, I, yt, yb, h, B, C, hf, Tt, Tb Z_top_min Z_bot_min]';
    section_properties_units = {'Mpa', 'Mpa', 'kN/m^3', '-', 'mm^3', 'mm^3', 'mm^2', 'mm^4', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm^3', 'mm^3'}';
    
    Tendon_properties_names = {'dp', 'e', 'Pi', 'Pe', 'nwire', 'fpu', 'db', 'n_tendon', 'Ap', 'A_tendon', 'Wp'}';
    Tendon_properties_values = [dp e P_i*1e-3 P_e*1e-3 n_w f_pu diam n_s An_ogd An_tendon W_p]';
    Tendon_properties_units = {'mm', 'mm', 'kN', 'kN', '-', 'Mpa', 'mm^2', '-', 'mm^2', 'mm^2', 'kg/m'}';
    
    Moment_kontrol_names = {'Mu', 'fi_Mn', 'Mcr', '1.2Mcr', 'Vu', 'fi_Vt'}';
    Moment_kontrol_values = [Mu fi_Mn Mcr Mcr_k V_u fi_Vt]';
    Moment_kontrol_units = {'kN.m', 'kN.m', 'kN.m', 'kN.m', 'kN', 'kN'}';
    
    Cost_name = {'Cost'};
    Cost_value = Cost;
    Cost_unit = {'$'};
    
    Rownames = [section_properties_names; Tendon_properties_names; Moment_kontrol_names; Cost_name];
    Constant = [section_properties_values; Tendon_properties_values; Moment_kontrol_values; Cost_value];
    Units = categorical([section_properties_units; Tendon_properties_units; Moment_kontrol_units; Cost_unit]);      
    Constant_section_properties = table(Constant, Units, 'RowNames',Rownames);
        
    %% Limit Zone For Tendon Centroid
    %**Midspan**%

    % e
    M_e = (-200:5:900)';
    M_size = size(M_e,1);

    % 1/Pi>fti
    M_1_bo_Pi_bu_fti(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fti(i,1) = (-1 + M_e(i) * yt / I * A) / (fti + Mo * 1e6 / Zt) / A; 
    end

    % 1/Pi>fci
    M_1_bo_Pi_bu_fci(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fci(i,1) = (1 + M_e(i) * yb / I * A) / (-fci + Mo * 1e6 / Zb) / A; 
    end

    % 1/Pi<fts
    M_1_bo_Pi_ku_fts(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + M_e(i) * yb / I * A) / (-fts + Mt * 1e6 / Zb) / A; 
    end

    % 1/Pi<fcs
    M_1_bo_Pi_ku_fcs(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + M_e(i) * yt / I * A) / (fcs + Mt * 1e6 / Zt) / A; 
    end

    % Pi
    M_Pi(M_size,4)=0;
    for i = 1:M_size
        M_Pi(i,1) = 1 / M_1_bo_Pi_bu_fti(i);
        M_Pi(i,2) = 1 / M_1_bo_Pi_bu_fci(i);
        M_Pi(i,3) = 1 / M_1_bo_Pi_ku_fts(i);
        M_Pi(i,4) = 1 / M_1_bo_Pi_ku_fcs(i);
    end

    M_Max = max([M_1_bo_Pi_bu_fti; M_1_bo_Pi_bu_fci; M_1_bo_Pi_ku_fts; M_1_bo_Pi_ku_fcs]);

    %***Graphic 1***%
    figure(1)
        plot(M_1_bo_Pi_bu_fti, M_e)
        title('Design of Beam with Constant Eccentricity Tendons - Midspan')
        ylabel('Eccentricity, mm')
        ylim([M_e(1) - 200 M_e(end) + 200])
        axis ij
        xlim([0 M_Max])
        grid on
        hold on
        plot(M_1_bo_Pi_bu_fci, M_e)
        plot(M_1_bo_Pi_ku_fts, M_e)
        plot(M_1_bo_Pi_ku_fcs, M_e)
        plot([0 M_Max], [e e])
        plot(1 / P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off
        
    %**Support**%
    % e
    S_e = (-200:5:900)';
    S_size = size(S_e,1);

    % 1/Pi>fti
    S_1_bo_Pi_bu_fti(S_size,1)=0;
    for i = 1:M_size
    S_1_bo_Pi_bu_fti(i,1) = (-1 + M_e(i) * yt / I * A) / (2 * fti * A);
    end

    % 1/Pi>fci
    S_1_bo_Pi_bu_fci(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_bu_fci(i,1) = (1 + M_e(i) * yb / I * A) / (-fci * A);
    end

    % 1/Pi<fts
    S_1_bo_Pi_ku_fts(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + M_e(i) * yb / I * A) / (-fts * A);
    end

    % 1/Pi<fts
    S_1_bo_Pi_ku_fcs(S_size,1)=0;
    for i = 1:S_size
        S_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + M_e(i) * yt / I * A) / (fcs * A);
    end

    % Pi
    S_Pi(S_size,4)=0;
    for i = 1:S_size
        S_Pi(i,1) = 1 / S_1_bo_Pi_bu_fti(i);
        S_Pi(i,2) = 1 / S_1_bo_Pi_bu_fci(i);
        S_Pi(i,3) = 1 / S_1_bo_Pi_ku_fts(i);
        S_Pi(i,4) = 1 / S_1_bo_Pi_ku_fcs(i);
    end

    S_Max = max([S_1_bo_Pi_bu_fti; S_1_bo_Pi_bu_fci; S_1_bo_Pi_ku_fts; S_1_bo_Pi_ku_fcs]);

    %***Graphic 2***%
    figure(2)
        plot(S_1_bo_Pi_bu_fti, S_e)
        title('Design of Beam with Constant Eccentricity Tendons - Support')
        ylabel('Eccentricity, mm')
        ylim([S_e(1) - 200 S_e(end) + 200])
        axis ij
        xlim([0 S_Max])
        grid on
        hold on
        plot(S_1_bo_Pi_bu_fci, S_e)
        plot(S_1_bo_Pi_ku_fts, S_e)
        plot(S_1_bo_Pi_ku_fcs, S_e)
        plot([0 S_Max], [e e])
        plot(1 / P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off

    %% Saving
    % Table
    if Opt_Method == 1
        save('ABC_Constant.mat','Constant_section_properties')
    elseif Opt_Method == 2
        save('BSO_Constant.mat','Constant_section_properties')
    else
        save('SSO_Constant.mat','Constant_section_properties')
    end

    % figure 1
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Constant Eccentricity Tendons - Midspan.png')
            delete('ABC Design of Beam with Constant Eccentricity Tendons - Midspan.png')
        end
        saveas(figure(1),'ABC Design of Beam with Constant Eccentricity Tendons - Midspan.png')
	elseif Opt_Method == 2
            if isfile('BSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
                delete('BSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'BSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
	else
            if isfile('SSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
                delete('SSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'SSO Design of Beam with Constant Eccentricity Tendons - Midspan.png')
    end

    % figure 2
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Constant Eccentricity Tendons - Support.png')
            delete('ABC Design of Beam with Constant Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'ABC Design of Beam with Constant Eccentricity Tendons - Support.png')
    elseif Opt_Method == 2
        if isfile('BSO Design of Beam with Constant Eccentricity Tendons - Support.png')
            delete('BSO Design of Beam with Constant Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'BSO Design of Beam with Constant Eccentricity Tendons - Support.png')
    else
        if isfile('SSO Design of Beam with Constant Eccentricity Tendons - Support.png')
            delete('SSO Design of Beam with Constant Eccentricity Tendons - Support.png')
        end
        saveas(figure(2),'SSO Design of Beam with Constant Eccentricity Tendons - Support.png')
    end

end

end

%% Design of Beam with Deflected Eccentricity Tendons
function [p_cost,penalty,fitness,sect] = DET(sect, write)
%% Required Input variables
global L w_l w_sdl gama_p beta1 fci fti fcs fts Section_Name fp_c fp_ci gama_conc An_ogd nu f_pu...
    Zt Zb A I yt yb h B C hf Tt Tb e dp n_w diam tend_pr W_p conc_pr n_s...
    Deflected_section_properties Opt_Method

%% Midspan moment for recomputed design loads
wo = (gama_conc * A * 1e-6); % (kN/m)
wd = w_sdl * (B + 2*C) * 1e-3; % (kN/m)
wl = w_l * (B + 2*C) * 1e-3; % (kN/m)

Mo = wo * L^2 / 8; % (kN.m)
Mop = wo * L^2 / 9; % (kN.m)
Md = (wd * L^2) / 8; % (kN.m)
Ml = (wl * L^2) / 8; % (kN.m)
Mt = Mo + Md + Ml; % (kN.m) %~ DDM
Mtp = (wo + wd + wl) * L^2 / 9; % (kN.m)
Mu = 1.4 * (Mo + Md) + 1.7 * Ml ; % (kN.m)

%% Minimum required values for the section modulus %~ DDM
Z_top_min= (-nu * Mop + Mo + Md + Ml) * 1e6 / (nu*fti - fcs); % (mm^3) %~ SDM
Z_bot_min= (-nu * Mop + Mo + Md + Ml) * 1e6 / (fts - nu*fci); % (mm^3) %~ SDM

if abs(Zt)>=abs(Z_top_min)
    pen_sec1 = 0;
 else
    pen_sec1 = abs(Z_top_min) / abs(Zt) - 1;
 end
if abs(Zb)>=abs(Z_bot_min)
	pen_sec2=0;
else
	pen_sec2 = abs(Z_bot_min) / abs(Zb) - 1;
end
pen_sec = pen_sec1 + pen_sec2;

%% Prestressing Force at Transfer Pi 
% fcci
fcci = fti - (yt / h) * (fti-fci); % (MPa)

% Pi
r2 = I / A; % (mm^2)
P_i = ((nu *(1 + e * yb / r2) / (- fts + Mt * 1e6 / Zb) / A) ^-1); % (N)

% emax
e_max = (fti - fcci) * (Zt / P_i) + Mop * 1e6 / P_i; % (mm)

% Pen_e
if abs(e) / abs(e_max) <= 1
    pen_e = 0;
else
    pen_e = abs(e) / abs(e_max) - 1;
end

%% Effective prestress force after all losses
P_e = nu * P_i ; % (N)

%% Specified yield strength of prestressing steel 
f_py = 0.85 * f_pu; % (MPa)

f_py_a = 0.82 * f_py; % (MPa)

f_pu_a = 0.74 * f_pu; % (MPa)

%%  Area of prestressing steel
A_tendon = P_i / min(f_py_a, f_pu_a); % (mm^2)

%% Number of strands
n_s = round(A_tendon / An_ogd + 0.5);
if n_s>20
    penn = n_s/20-1;
else
    penn=0;
end
%% New tendon
An_tendon = n_s * An_ogd; % (mm^2)
Pn_i = 0.82 * 0.85 * f_pu * An_tendon; % (N)
Pn_e = Pn_i * nu; % (N)

%% Flexural control
% Stress in prestressed reinforcement at nominal strength of component
ho_p = An_tendon / ((B + 2*C) * (dp)); % (birimsiz)
f_ps = f_pu * (1 - (gama_p / beta1) * (ho_p * f_pu) / fp_c); % (MPa)

% Depth of equivalent rectangular stress block
a = An_tendon * f_ps / (0.85 * fp_c * (B + 2*C)); % (mm)

% Nominal flexural strength at section
Mn = An_tendon * f_ps * (dp - a / 2) * 1e-6; % (kN.m)
fi_Mn = 0.9 * Mn; % (kN.m)

% pen_flexural
if abs(Mu) / abs(fi_Mn) <= 1
    pen_gocme = 0;
else
    pen_gocme = abs(Mu) / abs(fi_Mn) - 1;
end

%% Cracking control
Mcr = (0.35 * sqrt(fp_c) * Zb + Pn_e * (I / A / yb + e)) * 1e-6; % (kN.m)
Mcr_k = 1.2 * Mcr; % (kN.m)

if abs(Mcr_k) / abs(fi_Mn) <= 1
    pen_cr = 0;
else
    pen_cr = abs(Mcr_k) / abs(fi_Mn) - 1;
end

%% Shear Control
% Distance from extreme compression fiber to centroid of prestressing steel in tensionzone
d_p = max(0.8 * h, dp); % (mm)

% Web width
b_w = 12.5 * 25.4; % (mm)

% Shear forces
Vo = wo * L / 2; % (kN)
D_Vu = (wd + wl) * L / 2; % (kN)

% Moments
D_Mu = (wd + wl) * L^2 / 8; % (kN.m)
D_Mcr = Mcr - Mo; % (kN.m)

% Nominal shear strength provided by concrete when diagonal cracking results from combined shear and moment
V_ci1 = 0.05 * sqrt(fp_c) * b_w * d_p + Vo * 1e3 + D_Vu / D_Mu * D_Mcr * 1e3; % (N)
V_ci2 = 1.7 * sqrt(fp_c) * b_w * d_p; % (N)
V_ci = max(V_ci1, V_ci2) * 1e-3; % (kN)

% Nominal shear strength provided by concrete when diagonal cracking results from high principal tensile stress in the web 
V_cw = (0.29 * sqrt(fp_c) + 0.3 * P_e / A) * b_w * d_p * 1e-3; % (kN)

% Nominal shear strength provided by concrete
V_c = min(V_ci, V_cw); % (kN)

% Factored shear force at section
V_u = (1.4 * (wo + wd) + 1.7 * wl) * L / 2; % (kN)

% Area of shear reinforcement within spacing s=200
A_v1 = (V_u / 0.75 - V_c) * 1e3 * 200 / 420 / d_p; % (mm^2/200)
A_v2 = A_tendon / 80 * f_pu / 420 * 200 / d_p * sqrt(d_p / b_w); % (mm^2/200) 
A_v3 = 0.062 * sqrt(fp_c) * b_w * 200 / 420; % (mm^2/200)
A_v4 = 0.35 * b_w * 200 / 420; % (mm^2/200)
A_v = max([A_v1, A_v2, A_v3, A_v4]); % (mm^2/200)

% db_t = sqrt(A_v * 4 / pi) / 2; % (mm)

% Nominal shear strength provided by shear reinforcement
V_s = A_v / 200 * 420 * d_p * 1e-3; % (kN)

% pen_Shear
fi_Vt = 0.75 * (V_c + V_s); % (kN)

if abs(V_u) / abs(fi_Vt) <= 1
    pen_Shear = 0;
else
    pen_Shear = abs(V_u)/abs(fi_Vt)-1;
end

%% ACI Allowable Stresses(Midspan)
% ***Transmission***

% Pi/Ac
fb_MA_Pi = -P_i / A; % (MPa)
ft_MA_Pi = -P_i / A; % (MPa)

% Pie/Z
fb_MA_Pie = -P_i * e / Zb; % (MPa)
ft_MA_Pie = P_i * e / Zt; % (MPa)

% Mo/z
fb_MA_Mo = Mo * 1e6 / Zb; % (MPa)
ft_MA_Mo = -Mo * 1e6 / Zt; % (MPa)

% Total
fb_MA_t = fb_MA_Pi + fb_MA_Pie + fb_MA_Mo; % (MPa)
ft_MA_t = ft_MA_Pi + ft_MA_Pie + ft_MA_Mo; % (MPa)

% Allowable Sigma
% For fb
if fb_MA_t > 0
    fb_MA_se = fti; % (MPa)
else
    fb_MA_se = fci; % (MPa)
end
% For ft
if ft_MA_t > 0
    ft_MA_se = fti; % (MPa)
else
    ft_MA_se = fci; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_MA_se) >= abs(fb_MA_t)
    pen_MA1 = 0;
else
    pen_MA1 = abs(fb_MA_t) / abs(fb_MA_se) - 1;
end
% For ft
if abs(ft_MA_se) >= abs(ft_MA_t)
	pen_MA2 = 0;
else
	pen_MA2 = abs(ft_MA_t) / abs(ft_MA_se) - 1;
end
pen_MA = pen_MA1 + pen_MA2;
    
% ***Effective***

% Pe/Ac
fb_ME_Pe = -P_e / A; % (MPa)
ft_ME_Pe = -P_e / A; % (MPa)

% Pee/Z
fb_ME_Pee = -P_e * e / Zb; % (MPa)
ft_ME_Pee = P_e * e / Zt; % (MPa)

% Mt/z
fb_ME_Mt = Mt * 1e6 / Zb; % (MPa)
ft_ME_Mt = -Mt * 1e6 / Zt; % (MPa)

% Total
fb_ME_t = fb_ME_Pe + fb_ME_Pee + fb_ME_Mt; % (MPa)
ft_ME_t = ft_ME_Pe + ft_ME_Pee + ft_ME_Mt; % (MPa)

% Allowable Sigma
% For fb
if fb_ME_t > 0
    fb_ME_se = fts; % (MPa)
else
    fb_ME_se = fcs; % (MPa)
end
% for ft
if ft_ME_t > 0
    ft_ME_se = fts; % (MPa)
else
    ft_ME_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_ME_se) >= abs(fb_ME_t)
    pen_ME1=0;
else
    pen_ME1 = abs(fb_ME_t) / abs(fb_ME_se) - 1;
end
% For ft
if abs(ft_ME_se) >= abs(ft_ME_t)
	pen_ME2 = 0;
else
	pen_ME2 = abs(ft_ME_t) / abs(ft_ME_se) - 1;
end
pen_ME = pen_ME1 + pen_ME2;

%% Aci allowable stress(Draping)
% ***Transmission***

% Pi/Ac
fb_DA_Pi = -P_i / A; % (MPa)
ft_DA_Pi = -P_i / A; % (MPa)

% Pie/Z
fb_DA_Pie = -P_i * e / Zb; % (MPa)
ft_DA_Pie = P_i * e / Zt; % (MPa)

% Mt/z
fb_DA_Mo = Mop * 1e6 / Zb; % (MPa)
ft_DA_Mo = -Mop * 1e6 / Zt; % (MPa)

% Total
fb_DA_t = fb_DA_Pi + fb_DA_Pie + fb_DA_Mo; % (MPa) % Yeni Pi kullanýldý
ft_DA_t = ft_DA_Pi + ft_DA_Pie + ft_DA_Mo; % (MPa)

% Allowable Sigma
% for fb
if fb_DA_t > 0
    fb_DA_se = 2 * fti; % (MPa)
else
    fb_DA_se = fci; % (MPa)
end
% For ft
if ft_DA_t > 0
    ft_DA_se = 2 * fti; % (MPa)
else
    ft_DA_se = fci; % (MPa)
end

% pen_ACI_S
% For fb
if abs(fb_DA_se) >= abs(fb_DA_t)
    pen_DA1 = 0;
else
    pen_DA1 = abs(fb_DA_t) / abs(fb_DA_se) - 1;
end
% For ft
if abs(ft_DA_se) >= abs(ft_DA_t)
	pen_DA2 = 0;
else
	pen_DA2 = abs(ft_DA_t) / abs(ft_DA_se) - 1;
end
pen_DA = pen_DA1 + pen_DA2;

% ***Effective***

% Pe/Ac
fb_DE_Pe = -P_e / A; % (MPa)
ft_DE_Pe = -P_e / A; % (MPa)

% Pee/Z
fb_DE_Pee = -P_e * e / Zb; % (MPa)
ft_DE_Pee = P_e * e / Zt; % (MPa)

% Mt/z
fb_DE_Mt = Mtp * 1e6 / Zb; % (MPa)
ft_DE_Mt = -Mtp * 1e6 / Zt; % (MPa)

% Total
fb_DE_t = fb_DE_Pe + fb_DE_Pee + fb_DE_Mt; % (MPa)
ft_DE_t = ft_DE_Pe + ft_DE_Pee + ft_DE_Mt; % (MPa)

% Allowable Sigma
% For fb
if fb_DE_t > 0
    fb_DE_se = fts; % (MPa)
else
    fb_DE_se = fcs; % (MPa)
end
% For ft
if ft_DE_t > 0
    ft_DE_se = fts; % (MPa)
else
    ft_DE_se = fcs; % (MPa)
end

% pen_ACI_M
% For fb
if abs(fb_DE_se) >= abs(fb_DE_t)
    pen_DE1 = 0;
else
    pen_DE1 = abs(fb_DE_t) / abs(fb_DE_se) - 1;
end
% For ft
if abs(ft_DE_se) >= abs(ft_DE_t)
	pen_DE2 = 0;
else
	pen_DE2 = abs(ft_DE_t) / abs(ft_DE_se) - 1;
end
pen_DE = pen_DE1 + pen_DE2;

%% Costs 
% Cost of concrete
conc_cost = conc_pr * (A * 1e-6);

% Cost of tendon
tendon_cost = n_s * n_w * tend_pr * W_p;

% Cost of Reinforcement Steel
   % For a price of reinforcement steel 250 Euro/m^3, we have:
reinf_cost = 250 * A_v * 1e-6;

% Total Cost
Cost = conc_cost + tendon_cost + reinf_cost ; % Cost per metre

%% X Evaluate the Design

pen_str = pen_MA + pen_ME + pen_DA + pen_DE;
pen_UM = pen_gocme + pen_cr;
penalty = pen_sec + pen_e + pen_UM + pen_Shear + pen_str + penn;
p_cost = Cost * (1 + penalty)^2;
fitness = 1 / p_cost;

%% Displaying
if write == 1
    disp('---------------------------------------------------------')
    if Opt_Method == 1
        disp('<strong>********** ABC Design of Beam with Deflected Eccentricity Tendons **********</strong>')
    elseif Opt_Method == 2
        disp('<strong>********** BSO Design of Beam with Deflected Eccentricity Tendons **********</strong>')
    else
        disp('<strong>********** SSO Design of Beam with Deflected Eccentricity Tendons **********</strong>')
    end
    disp('<strong>Section Properties</strong>')
    fprintf('Selected Section Name : <strong>%s\n</strong>', Section_Name(:))
    fprintf('Zt = %7.4f mm^3\n', Zt)
    fprintf('Zb = %7.4f mm^3\n', Zb)
    fprintf('A  = %7.4f mm^2\n', A)
    fprintf('I  = %7.4f mm^4\n', I)
    fprintf('yt = %7.4f mm\n', yt)
    fprintf('yb = %7.4f mm\n', yb)
    fprintf('h  = %7.4f mm\n', h)
    fprintf('B  = %7.4f mm\n', B)
    fprintf('C  = %7.4f mm\n', C)
    fprintf('hf = %7.4f mm\n', hf)
    fprintf('Tt = %7.4f mm\n', Tt)
    fprintf('Tb = %7.4f mm\n', Tb)
    fprintf('Position of prestress tendon = %7.4f mm\n', dp)
    fprintf('                Eccentricity = %7.4f mm\n', e)
    fprintf('  Initial Prestressing Force = %7.4f kN\n', P_i*1e-3)
    fprintf('  Effective prestress force  = %7.4f kN\n', P_e*1e-3)
    
    disp('---------------------------------------------------------')
    disp('<strong>Tendon Properties</strong>')
    fprintf('  nwire = %7.4f \n', n_w)
    fprintf('    fpu = %7.4f MPa\n', f_pu)
    fprintf('     db = %7.4f mm^2\n', diam)
    fprintf('ntendon = %7.4f \n', n_s)
    fprintf('     Ap = %7.4f mm^2\n', An_ogd)
    fprintf('Atendon = %7.4f mm^2\n', An_tendon)
    fprintf('     Wp = %7.4f kg/m\n', W_p)

    disp('---------------------------------------------------------')
    disp('  <strong>I. Minimum required values for the section modulus</strong>')
    if Zt >= abs(Z_top_min) && Zb >= abs(Z_bot_min)
        fprintf( 'Zt = %3.4g mm > Z_top_min = %3.4g mm\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm > Z_bot_min = %3.4g mm\n', Zb, abs(Z_bot_min))
        disp('<strong>OK</strong>');
    else    
        fprintf( 'Zt = %3.4g mm >= Z_top_min = %3.4g mm olmasi lazim\n', Zt, abs(Z_top_min))
        fprintf( 'Zb = %3.4g mm >= Z_bot_min = %3.4g mm olmasi lazim\n', Zb, abs(Z_bot_min))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>II. Flexural Control</strong>')
    if Mu / fi_Mn <= 1
        fprintf('Mu = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mu = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mu, fi_Mn)
        disp('<strong>NOT OK</strong>')
    end

    disp(' ')
    disp('  <strong>III. Cracking Control</strong>')
    if Mcr_k / fi_Mn <= 1
        fprintf('Mcr_k = %3.4g kN.m <= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>OK</strong>');
    else
        fprintf('Mcr_k = %3.4g kN.m >= fi_Mn = %3.4g kN.m\n', Mcr_k, fi_Mn)
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>IV. Shear Control</strong>')
    if  V_u / fi_Vt <= 1
        fprintf('V_u = %3.4g kN <= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>OK</strong>');
    else
        fprintf('V_u = %3.4g kN >= fi_Vt = %3.4g kN\n', V_u, fi_Vt)
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('  <strong>V. Aci Allowable Stresses</strong>')
    disp(' <strong>V.1 Midspan</strong>')
    disp('<strong>V.1.1 Transmission</strong>')

    % For fb
    if  abs(fb_MA_t) / abs(fb_MA_se) <= 1
        fprintf('abs(fb_MA_t) = %3.4g MPa <= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_MA_t) = %3.4g MPa >= abs(fb_MA_se) = %3.4g MPa\n', abs(fb_MA_t), abs(fb_MA_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if  abs(ft_MA_t) / abs(ft_MA_se) <= 1
        fprintf('abs(ft_MA_t) = %3.4g MPa <= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_MA_t) = %3.4g MPa >= abs(ft_MA_se) = %3.4g MPa\n', abs(ft_MA_t), abs(ft_MA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
     disp('<strong>V.1.2 Effective</strong>')
    % For fb
    if  abs(fb_ME_t) / abs(fb_ME_se) <= 1
        fprintf('abs(fb_ME_t) = %3.4g MPa <= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_ME_t) = %3.4g MPa >= abs(fb_ME_se) = %3.4g MPa\n', abs(fb_ME_t), abs(fb_ME_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_ME_t) / abs(ft_ME_se) <= 1
        fprintf('abs(ft_ME_t) = %3.4g MPa <= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_ME_t) = %3.4g MPa >= abs(ft_ME_se) = %3.4g MPa\n', abs(ft_ME_t), abs(ft_ME_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp(' <strong>V.2 Draping</strong>')
    disp('<strong>V.2.1 Transmission</strong>')
    % For fb
    if  abs(fb_DA_t) / abs(fb_DA_se) <= 1
        fprintf('abs(fb_DA_t) = %3.4g MPa <= abs(fb_DA_se) = %3.4g MPa\n', abs(fb_DA_t), abs(fb_DA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_DA_t) = %3.4g MPa >= abs(fb_DA_se) = %3.4g MPa\n', abs(fb_DA_t), abs(fb_DA_se))
        disp('<strong>NOT OK</strong>');
    end
    % ft icin
    if  abs(ft_DA_t) / abs(ft_DA_se) <= 1
        fprintf('abs(ft_DA_t) = %3.4g MPa <= abs(ft_DA_se) = %3.4g MPa\n', abs(ft_DA_t), abs(ft_DA_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_DA_t) = %3.4g MPa >= abs(ft_DA_se) = %3.4g MPa\n', abs(ft_DA_t), abs(ft_DA_se))
        disp('<strong>NOT OK</strong>');
    end

    disp(' ')
    disp('<strong>V.2.2 Effective</strong>')
    % For fb
    if  abs(fb_DE_t) / abs(fb_DE_se) <= 1
        fprintf('abs(fb_DE_t) = %3.4g MPa <= abs(fb_DE_se) = %3.4g MPa\n', abs(fb_DE_t), abs(fb_DE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(fb_DE_t) = %3.4g MPa >= abs(fb_DE_se) = %3.4g MPa\n', abs(fb_DE_t), abs(fb_DE_se))
        disp('<strong>NOT OK</strong>');
    end
    % For ft
    if abs(ft_DE_t) / abs(ft_DE_se) <= 1
        fprintf('abs(ft_DE_t) = %3.4g MPa <= abs(ft_DE_se) = %3.4g MPa\n', abs(ft_DE_t), abs(ft_DE_se))
        disp('<strong>OK</strong>');
    else
        fprintf('abs(ft_DE_t) = %3.4g MPa >= abs(ft_DE_se) = %3.4g MPa\n', abs(ft_DE_t), abs(ft_DE_se))
        disp('<strong>NOT OK</strong>');
    end

    disp('---------------------------------------------------------');
    disp('<strong>Penalty</strong>')
    fprintf('      Penalty from Section = %7.4f \n',pen_sec)
    fprintf(' Penalty from Eccentricity = %7.4f \n', pen_e)
    fprintf(' Penalty from Moment_Check = %7.4f \n', pen_UM)
    fprintf('  Penalty from Shear_Check = %7.4f \n', pen_Shear)
    fprintf('Penalty from Tendon stress = %7.4f \n', pen_str)

    %% Table
    section_properties_names = {'fcp', 'fcpi', 'Gama conc', 'Ratio nu', 'Zt', 'Zb', 'A', 'I', 'yt', 'yb', 'h', 'B', 'C', 'hf', 'Tt', 'Tb','Zt_min', 'Zb_min'}';
    section_properties_values = [fp_c, fp_ci, gama_conc, nu, Zt, Zb, A, I, yt, yb, h, B, C, hf, Tt, Tb Z_top_min Z_bot_min]';
    section_properties_units = {'Mpa', 'Mpa', 'kN/m^3', '-', 'mm^3', 'mm^3', 'mm^2', 'mm^4', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm', 'mm^3', 'mm^3'}';
    
    Tendon_properties_names = {'dp', 'e', 'Pi', 'Pe', 'nwire', 'fpu', 'db', 'n_tendon', 'Ap', 'A_tendon', 'Wp'}';
    Tendon_properties_values = [dp e P_i*1e-3 P_e*1e-3 n_w f_pu diam n_s An_ogd An_tendon W_p]';
    Tendon_properties_units = {'mm', 'mm', 'kN', 'kN', '-', 'Mpa', 'mm^2', '-', 'mm^2', 'mm^2', 'kg/m'}';
    
    Moment_kontrol_names = {'Mu', 'fi_Mn', 'Mcr', '1.2Mcr', 'Vu', 'fi_Vt'}';
    Moment_kontrol_values = [Mu fi_Mn Mcr Mcr_k V_u fi_Vt]';
    Moment_kontrol_units = {'kN.m', 'kN.m', 'kN.m', 'kN.m', 'kN', 'kN'}';
    
    Cost_name = {'Cost'};
    Cost_value = Cost;
    Cost_unit = {'$'};
    
    Rownames = [section_properties_names; Tendon_properties_names; Moment_kontrol_names; Cost_name];
    Deflected = [section_properties_values; Tendon_properties_values; Moment_kontrol_values; Cost_value];
    Units = categorical([section_properties_units; Tendon_properties_units; Moment_kontrol_units; Cost_unit]);   
    Deflected_section_properties = table(Deflected, Units, 'RowNames',Rownames);
     
    %% Limit Zone For Tendon Centroid
    %**Midspan**%

    % e
    M_e = (-200:5:900)';
    M_size = size(M_e,1);

    % 1/Pi>fti
    M_1_bo_Pi_bu_fti(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fti(i,1) = (-1 + M_e(i) * yt / I * A) / (fti + Mo * 1e6 / Zt) / A; 
    end

    % 1/Pi>fci
    M_1_bo_Pi_bu_fci(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_bu_fci(i,1) = (1 + M_e(i) * yb / I * A) / (-fci + Mo * 1e6 / Zb) / A; 
    end

    % 1/Pi<fts
    M_1_bo_Pi_ku_fts(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + M_e(i) * yb / I * A) / (-fts + Mt * 1e6 / Zb) / A; 
    end

    % 1/Pi<fcs
    M_1_bo_Pi_ku_fcs(M_size,1)=0;
    for i = 1:M_size
        M_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + M_e(i) * yt / I * A) / (fcs + Mt * 1e6 / Zt) / A; 
    end

    % Pi
    M_Pi(M_size,4)=0;
    for i = 1:M_size
        M_Pi(i,1) = 1 / M_1_bo_Pi_bu_fti(i);
        M_Pi(i,2) = 1 / M_1_bo_Pi_bu_fci(i);
        M_Pi(i,3) = 1 / M_1_bo_Pi_ku_fts(i);
        M_Pi(i,4) = 1 / M_1_bo_Pi_ku_fcs(i);
    end

    M_Max = max([M_1_bo_Pi_bu_fti; M_1_bo_Pi_bu_fci; M_1_bo_Pi_ku_fts; M_1_bo_Pi_ku_fcs]);

    %***Graphic***%
    figure(1)
        plot(M_1_bo_Pi_bu_fti, M_e)
        title('Design of Beam with Deflected Eccentricity Tendons - Midspan')
        ylabel('Eksentrisite, mm')
        ylim([M_e(1) - 200 M_e(end) + 200])
        axis ij
        xlim([0 M_Max])
        grid on
        hold on
        plot(M_1_bo_Pi_bu_fci, M_e)
        plot(M_1_bo_Pi_ku_fts, M_e)
        plot(M_1_bo_Pi_ku_fcs, M_e)
        plot([0 M_Max], [e e])
        plot(1 / P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off
        
    %**Draping**%
    % e
    D_e = (-200:5:900)';
    D_size = size(D_e,1);

    % 1/Pi>fti
    D_1_bo_Pi_bu_fti(M_size,1)=0;
    for i = 1:D_size
        D_1_bo_Pi_bu_fti(i,1) = (-1 + D_e(i) * yt / I * A) / (fti + Mop * 1e6 / Zt) / A; %~ SDM
    end

    % 1/Pi>fci
    D_1_bo_Pi_bu_fci(D_size,1)=0;
    for i = 1:D_size
        D_1_bo_Pi_bu_fci(i,1) = (1 + D_e(i) * yb / I * A) / (-fci + Mop * 1e6 / Zb) / A; %~ SDM
    end

    % 1/Pi<fts
    D_1_bo_Pi_ku_fts(D_size,1)=0;
    for i = 1:D_size
        D_1_bo_Pi_ku_fts(i,1) = 0.85 * (1 + D_e(i) * yb / I * A) / (-fts + Mtp * 1e6 / Zb) / A; %~ SDM
    end

    % 1/Pi<fts
    D_1_bo_Pi_ku_fcs(D_size,1)=0;
    for i = 1:D_size
        D_1_bo_Pi_ku_fcs(i,1) = 0.85 * (-1 + D_e(i) * yt / I * A) / (fcs + Mtp * 1e6 / Zt) / A; %~ SDM
    end

    % Pi
    D_Pi(D_size,4)=0;
    for i = 1:D_size
        D_Pi(i,1) = 1 / D_1_bo_Pi_bu_fti(i);
        D_Pi(i,2) = 1 / D_1_bo_Pi_bu_fci(i);
        D_Pi(i,3) = 1 / D_1_bo_Pi_ku_fts(i);
        D_Pi(i,4) = 1 / D_1_bo_Pi_ku_fcs(i);
    end

    D_Max = max([D_1_bo_Pi_bu_fti; D_1_bo_Pi_bu_fci; D_1_bo_Pi_ku_fts; D_1_bo_Pi_ku_fcs]);

    %***Graphic***%
    figure(2)
        plot(D_1_bo_Pi_bu_fti, D_e)
        title('Design of Beam with Deflected Eccentricity Tendons - Draping')
        ylabel('Eksentrisite, mm')
        ylim([D_e(1) - 200 D_e(end) + 200])
        axis ij
        xlim([0 D_Max])
        grid on
        hold on
        plot(D_1_bo_Pi_bu_fci, D_e)
        plot(D_1_bo_Pi_ku_fts, D_e)
        plot(D_1_bo_Pi_ku_fcs, D_e)
        plot([0 D_Max], [e e])
        plot(1 / P_i, e, '.', 'color', 'r', 'MarkerSize',20)
        legend('fti', 'fci', 'fts', 'fcs', 'Max', 'PS')
        hold off
        
    %% Saving
    % Table
    if Opt_Method == 1
        save('ABC_Deflected.mat','Deflected_section_properties')
    elseif Opt_Method == 2
        save('BSO_Deflected.mat','Deflected_section_properties')
    else
        save('SSO_Deflected.mat','Deflected_section_properties')
    end

    % figure 1
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
            delete('ABC Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
        end
        saveas(figure(1),'ABC Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
	elseif Opt_Method == 2
            if isfile('BSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
                delete('BSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'BSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
	else
            if isfile('SSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
                delete('SSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
            end
        saveas(figure(1),'SSO Design of Beam with Deflected Eccentricity Tendons - Midspan.png')
    end

    % figure 2
    if Opt_Method == 1
        if isfile('ABC Design of Beam with Deflected Eccentricity Tendons - Draping.png')
            delete('ABC Design of Beam with Deflected Eccentricity Tendons - Draping.png')
        end
        saveas(figure(2),'ABC Design of Beam with Deflected Eccentricity Tendons - Draping.png')
    elseif Opt_Method == 2
        if isfile('BSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
            delete('BSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
        end
        saveas(figure(2),'BSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
    else
        if isfile('SSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
            delete('SSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
        end
        saveas(figure(2),'SSO Design of Beam with Deflected Eccentricity Tendons - Draping.png')
    end
end
end
