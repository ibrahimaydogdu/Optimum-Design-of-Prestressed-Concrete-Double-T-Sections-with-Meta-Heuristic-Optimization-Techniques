%% Clean screen and workspace
clear;
clc;

%% Table for Variable Eccentricity Tendons Beam design
load ABC_Variable.mat
var_ABC = Variable_section_properties;
var_ABC.Properties.VariableNames{'Variable'} = 'ABC';
clear ABC_Variable.mat Variable_section_properties

load BSO_Variable.mat
var_BSO = Variable_section_properties;
var_BSO.Properties.VariableNames{'Variable'} = 'BSO';
clear BSO_Variable.mat Variable_section_properties

load SSO_Variable.mat
var_SSO = Variable_section_properties;
var_SSO.Properties.VariableNames{'Variable'} = 'SSO';
clear SSO_Variable.mat Variable_section_properties

% Merging
var_table = [var_ABC(:,1),var_BSO(:,1),var_SSO];
var_table.Properties.RowNames
writetable(var_table,'summary.xlsx','Sheet',1,'Range','A2')

%% Table for Constant Eccentricity Tendons Beam design
load ABC_Constant.mat
const_ABC = Constant_section_properties;
const_ABC.Properties.VariableNames{'Constant'} = 'ABC';
clear ABC_Constant.mat Constant_section_properties

load BSO_Constant.mat
const_BSO = Constant_section_properties;
const_BSO.Properties.VariableNames{'Constant'} = 'BSO';
clear BSO_Constant.mat Constant_section_properties

load SSO_Constant.mat
const_SSO = Constant_section_properties;
const_SSO.Properties.VariableNames{'Constant'} = 'SSO';
clear SSO_Constant.mat Constant_section_properties

% Merging
const_table = [const_ABC(:,1),const_BSO(:,1),const_SSO];
writetable(const_table,'summary.xlsx','Sheet',2,'Range','A2')

%% Table for Deflected Eccentricity Tendons Beam design
load ABC_Deflected.mat
deflec_ABC = Deflected_section_properties;
deflec_ABC.Properties.VariableNames{'Deflected'} = 'ABC';
clear ABC_Deflected.mat Deflected_section_properties

load BSO_Deflected.mat
deflec_BSO = Deflected_section_properties;
deflec_BSO.Properties.VariableNames{'Deflected'} = 'BSO';
clear BSO_Deflected.mat Deflected_section_properties

load SSO_Deflected.mat
deflec_SSO = Deflected_section_properties;
deflec_SSO.Properties.VariableNames{'Deflected'} = 'SSO';
clear SSO_Deflected.mat Deflected_section_properties

% Merging
deflec_table = [deflec_ABC(:,1),deflec_BSO(:,1),deflec_SSO];
writetable(deflec_table,'summary.xlsx','Sheet',3,'Range','A2')

%% For Variable
load ABC_Variable_history.mat
load BSO_Variable_history.mat
load SSO_Variable_history.mat

% Grafic
figure(1)

% For Abc
plot(AVx, AVy,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'MarkerSize',5)

title('Search History of the optimization Algorithm - Variable')
xlabel('Iteration')
ylabel('Cost($)')
grid on
hold on

% For BSO
plot(BVx, BVy,'-r^','LineWidth',2,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',5)

% For SSO
plot(SVx, SVy,'-go','LineWidth',2,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g',...
    'MarkerSize',5)

hold off
legend('ABC', 'BSO', 'SSO')

% Save
saveas(figure(1),'Summary_Variable.png')

%% For COnstant
load ABC_Constant_history.mat
load BSO_Constant_history.mat
load SSO_Constant_history.mat

% Grafic
figure(2)

% For Abc
plot(ACx, ACy,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'MarkerSize',5)

title('Search History of the optimization Algorithm - Constant')
xlabel('Iteration')
ylabel('Cost($)')
grid on
hold on

% For BSO
plot(BCx, BCy,'-r^','LineWidth',2,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',5)

% For SSO
plot(SCx, SCy,'-go','LineWidth',2,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g',...
    'MarkerSize',5)

hold off
legend('ABC', 'BSO', 'SSO')

% Save
saveas(figure(2),'Summary_Constant.png')

%% For Defleted
load ABC_Deflected_history.mat
load BSO_Deflected_history.mat
load SSO_Deflected_history.mat

% Grafic
figure(3)
% For Abc
plot(ADx, ADy,'-bs','LineWidth',2,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor','b',...
    'MarkerSize',5)

title('Search History of the optimization Algorithm - Defleted')
xlabel('Iteration')
ylabel('Cost($)')
grid on
hold on

% For BSO
plot(BDx, BDy,'-r^','LineWidth',2,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',5)

% For SSO
plot(SDx, SDy,'-go','LineWidth',2,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor','g',...
    'MarkerSize',5)

hold off
legend('ABC', 'BSO', 'SSO')

% Save
saveas(figure(3),'Summary_Defleted.png')