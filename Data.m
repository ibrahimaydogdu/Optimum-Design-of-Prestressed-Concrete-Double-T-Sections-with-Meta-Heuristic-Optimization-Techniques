%% Global variable
global units Section_name ub lb n conc_grad nu_M section tendon dpten 

%% ***** Type or grade of concrete (ffc, ffci) *****
% *** SI units ***
ffc_SI = (15:5:80)'; % (MPa)
ffci_SI = 0.8*ffc_SI; % (MPa)

% *** US units ***
ffc_US = ffc_SI / 6.895e-3; % (psi)
ffci_US = ffci_SI / 6.895e-3; % (psi)

% Concrete price
conc_pr=[170 176 184 192 206 220 250 271 275 280 295 307 320 335]';

% *** Concrete grade ***
conc_grad_SI_M = [ffc_SI ffci_SI conc_pr];
conc_grad_US_M = [ffc_US ffci_US conc_pr];

%% Self weight of the concrete (Sabit olmasi lazim mi?)
% global g_conc_SI_M g_conc_US_M
% % *** SI units ***
% g_conc_SI_M = (15:0.5:30)'; % (kN/m^3)
% 
% % *** US units ***
% g_conc_US_M = g_conc_SI_M / 16.03; % (pcf)

%% ***** Ratio nu *****
nu_M = (0.75:0.01:0.85)';

%% Prestressing Tendon Type and Area
% Tendon(i)=[n_w,fpu,diam_p,Ap,Wp]
Tendon_US_M(1,1:5) = [7, 270, 3/8, 0.085, 0.29];
Tendon_US_M(2, 1:5) = [7, 270, 7/16, 0.115, 0.40];
Tendon_US_M(3, 1:5) = [7, 270, 1/2, 0.153, 0.52];
Tendon_US_M(4, 1:5) = [7, 270, 1/2, 0.167, 0.53]; % The 1/2 in. special strand has a larger actual diameter than the 1/2 in. regular strand.
Tendon_US_M(5, 1:5) = [7, 270, 9/16, 0.192, 0.65];
Tendon_US_M(6, 1:5) = [7, 270, 3/5, 0.217, 0.74];
Tendon_US_M(7, 1:5) = [7, 250, 1/4, 0.036, 0.12];
Tendon_US_M(8, 1:5) = [7, 250, 5/16, 0.058, 0.20];
Tendon_US_M(9, 1:5) = [7, 250, 3/8, 0.080, 0.27];
Tendon_US_M(10, 1:5) = [7, 250, 7/16, 0.108, 0.37];
Tendon_US_M(11, 1:5) = [7, 250, 1/2, 0.144, 0.49];
Tendon_US_M(12, 1:5) = [7, 250, 3/5, 0.216, 0.74];
Tendon_US_M(13, 1:5) = [3, 250, 1/4, 0.036, 0.12];
Tendon_US_M(14, 1:5) = [3, 250, 5/16, 0.058, 0.20];
Tendon_US_M(15, 1:5) = [3, 250, 3/8, 0.075, 0.26];
Tendon_US_M(16, 1:5) = [4, 250, 7/16, 0.106, 0.36];
Tendon_US_M(17, 1:5) = [1, 279, 0.105, 0.0087, 0.030];
Tendon_US_M(18, 1:5) = [1, 273, 0.120, 0.0114, 0.039];
Tendon_US_M(19, 1:5) = [1, 268, 0.135, 0.0143, 0.049];
Tendon_US_M(20, 1:5) = [1, 263, 0.148, 0.0173, 0.059];
Tendon_US_M(21, 1:5) = [1, 259, 0.162, 0.0206, 0.070];
Tendon_US_M(22, 1:5) = [1, 255, 0.177, 0.0246, 0.083];
Tendon_US_M(23, 1:5) = [1, 250, 0.192, 0.0289, 0.098];
Tendon_US_M(24, 1:5) = [1, 250, 0.196, 0.0302, 0.10];
Tendon_US_M(25, 1:5) = [1, 240, 0.250, 0.0491, 0.17];
Tendon_US_M(26, 1:5) = [1, 235, 0.276, 0.0598, 0.20];
Tendon_US_M(27, 1:5) = [1, 145, 3/8, 0.442, 1.50];
Tendon_US_M(28, 1:5) = [1, 145, 7/8, 0.601, 2.04];
Tendon_US_M(29, 1:5) = [1, 145, 1, 0.785, 2.67];
Tendon_US_M(30, 1:5) = [1, 145, 1+1/8, 0.994, 3.38];
Tendon_US_M(31, 1:5) = [1, 145, 1+1/4, 1.227, 4.17];
Tendon_US_M(32, 1:5) = [1, 145, 1+3/8, 1.485, 5.05];
Tendon_US_M(33, 1:5) = [1, 160, 3/8, 0.442, 1.50];
Tendon_US_M(34, 1:5) = [1, 160, 7/8, 0.601, 2.04];
Tendon_US_M(35, 1:5) = [1, 160, 1, 0.785, 2.67];
Tendon_US_M(36, 1:5) = [1, 160, 1+1/8, 0.994, 3.38];
Tendon_US_M(37, 1:5) = [1, 160, 1+1/4, 1.227, 4.17];
Tendon_US_M(38, 1:5) = [1, 160, 1+3/8, 1.485, 5.05];
Tendon_US_M(39, 1:5) = [1, 157, 5/8, 0.28, 0.98];
Tendon_US_M(40, 1:5) = [1, 150, 1, 0.85, 3.01];
Tendon_US_M(41, 1:5) = [1, 160, 1, 0.85, 3.01]; % Verify availability before specifying
Tendon_US_M(42, 1:5) = [1, 150, 1+1/4, 1.25, 4.39];
Tendon_US_M(43, 1:5) = [1, 160, 1+1/4, 1.25, 4.39]; % Verify availability before specifying
Tendon_US_M(44, 1:5) = [1, 150, 1+3/8, 1.58, 5.56];
Tendon_US_M(:, 6) = [720;720;720;720;720;720;680;680;680;780;780;780;750;750;750;750;750;750;720;720;750;750;750;750;735;800;600;800;800;800;820;1020;1030;1030;1020;1800;1800;1800;1800;1800;1800;1800;1800;1800];

% *** Sort by weight ***
Tendon_US_M = sortrows(Tendon_US_M, [1 5]);
indx = Tendon_US_M(:,1) < 6;
Tendon_US_M(indx,:) = [];
[nbr_rows_tend, nbr_clms_tend] = size(Tendon_US_M);
[Tendon_US_M_sorted, IX] = sortrows(Tendon_US_M, 5);

% *** SI units ***
conv_tend_US_to_SI = [1, 6.895, 25.4, 645.2, 1.488, 1e-3];
Tendon_SI_M_sorted = zeros(nbr_rows_tend, nbr_clms_tend);
for i = 1:nbr_rows_tend
    Tendon_SI_M_sorted(i, :) = Tendon_US_M_sorted(i, :) .*conv_tend_US_to_SI;
end

%% Position of prestress tendon

% *** SI units ***
dp_US_M = (0.1:0.01:0.9)';

% *** US units ***
dp_SI_M = dp_US_M;

%% ***** PCI Standart Double T Section *****

% PCI Design Handbook, 6th ed. 
% Birimler icin in ve mm
% A      : Kesit alani 
% I      : Atalet moment
% yt     : Ustten agirlik merkezine mesafe
% yb     : Alttan agirlik merkezine mesafe 
% Zt     : section modulus of top fibre
% Zb     : section modulus of bottom fibre
% sec(i)=[Zt,Zb,A,I,yt,yb,h,B,C,hf,Tt,Tb] 
sec1 = [1001, 315, 287, 2872, 2.87, 9.13, 12, 48, 24, 2, 5.75, 3.75];  % 8DT12   untoped

sec2 = [4508/3.49,4508/10.51,306,4508,3.49,10.51,14,48,24,2,5.75,3.75];  % 8DT14   untoped

sec3 = [6634/4.07 ,6634/11.93 , 325, 6634, 4.07, 11.93, 16, 48, 24, 2, 5.75, 3.75];  % 8DT16   untoped

sec4 = [9300/4.73, 9300/13.27, 344, 9300, 4.73, 13.27, 18, 48, 24, 2, 5.75, 3.75];  % 8DT18   untoped

sec5 = [12551/5.41, 12551/14.59, 363, 12551, 5.41, 14.59, 18, 48, 24, 2, 5.75, 3.75];  % 8DT20   untoped

sec6 = [20985/6.85, 20985/17.15, 401, 20985, 6.85, 17.15, 24, 48, 24, 2, 5.75, 3.75];  % 8DT24   untoped

sec7 = [5140, 2615, 567, 55464, 10.79, 21.21, 32, 48, 24, 2, 7.75, 4.75];  % 8DT32   untoped

sec8 = [3607, 1264, 449, 22469, 6.23, 17.77, 24, 60, 30, 2, 5.75, 3.75];  % 10DT24   untoped

sec9 = [5960, 2717, 615, 59720, 10.02, 21.98, 32, 60, 30, 2, 7.75, 4.75];  % 10DT32   untoped

sec10 = [5772, 2205, 640, 44563, 7.79, 20.21, 28, 72, 36, 2, 7.75, 4.75];  % 12DT28   untoped

sec11 = [6986, 2840, 690, 64620, 9.25, 22.75, 32, 72, 36, 2, 7.75, 4.75];  % 12DT32   untoped

sec12 = [5379, 1514, 689, 30716, 5.71, 20.29, 26, 60, 30, 4, 5.75, 3.75]; % 10DT26  Pretopped

sec13 = [9046, 3222, 855, 80780, 8.93, 25.07, 34, 60, 30, 4, 7.75, 4.75]; % 10DT34  Pretopped

sec14 = [8498, 2615, 928, 59997, 7.06, 22.94, 30, 72, 36, 4, 7.75, 4.75]; % 12DT30  Pretopped

sec15 = [10458, 3340, 978, 86072, 8.23, 25.77, 34, 72, 36, 4, 7.75, 4.75]; % 12DT34  Pretopped

sec16 = [8618, 2688, 1078, 53280, 6.18, 19.82, 26, 90, 45, 4, 9, 7.25]; % 15DT26  Pretopped

sec17 = [10838, 3457, 1133, 78625, 7.25, 22.75, 30, 90, 45, 4, 9, 6.875]; % 15DT30  Pretopped

sec18 = [13121,4274,1185,109621,8.32,25.65,34,90,45,4,9,6.50]; % 15DT34  Pretopped

Section_Name = ['8DT12 '; '8DT14 '; '8DT16 '; '8DT18 '; '8DT20 '; '8DT24 '; '8DT32 ';...
    '10DT24'; '10DT32'; '12DT28'; '12DT32'; '10DT26'; '10DT34'; '12DT30'; '12DT34';...
    '15DT26'; '15DT30'; '15DT34'];

% *** Selection of suitable section ***
% Matrix of all standart double T section, Each row represents properties of a double T section.
sec_all = [sec1; sec2; sec3; sec4; sec5; sec6; sec7; sec8; sec9; sec10; sec11;...
    sec12; sec13; sec14; sec15; sec16; sec17; sec18];

% sorts the elements of sort_sec matrix in ascending order along the A (Area of the section) array
[sort_sec_US_M, ind] = sortrows(sec_all,3);

%  Matrix row and column number
[nbr_rows_sec, nbr_clms_sec] = size(sort_sec_US_M);

Section_Name_Sorted=strings(nbr_rows_sec,1);
% Section_Name_Sorted(nbr_rows_sec, 1) = 0;
 for i=1:nbr_rows_sec
     Section_Name_Sorted(i)=Section_Name(ind(i),:);
 end
 
% *** Conversion Section properties from US to SI units *** %

% conversion matrix for Zt Zb A I yt yb h B C hf Tt Tb
conv_sec_US_to_SI = [16387, 16387, 645.2, 416231, 25.4, 25.4, 25.4, 25.4, 25.4, 25.4, 25.4, 25.4];

% conversion of Zt Zb A I yt yb h B C hf Tt Tb
sort_sec_SI_M = zeros(nbr_rows_sec, nbr_clms_sec);
for i = 1:nbr_rows_sec
    sort_sec_SI_M(i, :) = sort_sec_US_M(i, :) .* conv_sec_US_to_SI; 
end

%% Design variables of the optimization problem

%number of design varibles
nofdesignv = 4;
% All design variables are considered as discrete
% Define lower and upper boundries according to the possibilities
ub = zeros(nofdesignv, 1); lb = zeros(nofdesignv, 1);

% 2.2 Ratio of instantaneous prestress losses nu (11 possibilities)
ub(1) = numel(nu_M); lb(1) = 1; 

% 2.3 Double T beam sections (18 possibilities)
ub(2) = nbr_rows_sec; lb(2) = 1;       

% 2.4 Prestressing Tendon Type (41 possibilities)
ub(3) = nbr_rows_tend; lb(3) = 1;

% 2.5 Position of prestress tendon (701 possibilities)
ub(4) = numel(dp_US_M); lb(4) = 1;

% 2.5 Type or grade of concrete (14 possibilities)
ub(5) = numel(ffc_US); lb(5) = 1; 

% 2.6 Self weight of the concrete. (kN/m^3) (31 possibilities)
% ub(6) = numel(g_conc_US_M); lb(6) = 1; % Self weight of the concrete is taking constant 25kN/m^3     


%% Unit choice

if units == 1
    conc_grad = conc_grad_US_M;
    section = sort_sec_US_M;
    tendon = Tendon_US_M_sorted;
    dpten = dp_US_M;
else
    conc_grad = conc_grad_SI_M;
    section = sort_sec_SI_M;
    tendon = Tendon_SI_M_sorted;
    dpten = dp_SI_M;
end
Section_name = Section_Name_Sorted;
n = nofdesignv;

