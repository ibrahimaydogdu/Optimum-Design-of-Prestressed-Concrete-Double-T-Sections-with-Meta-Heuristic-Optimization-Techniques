%% 1. Clean
clear; close all; clc;

%% 2. Starting
global units Opt_Method design_Method write

disp ('Optimization program for design of prestressed double T beam');

units = input('Enter units system( US = 1, SI = 2): ');

% Units control
if units ~= 1 && units ~= 2
    disp('Please check units system')
    return
end

Opt_Method = input('Enter optimization method(ABC=1, BSO=2, SSO=3): ');

% Opt_Method control
if Opt_Method ~= 1 && Opt_Method ~= 2 && Opt_Method ~= 3
    disp('Please check optimization method')
    return
end

design_Method = input('Enter design method(Variable... = 1, Constant...= 2, Deflected...= 3): ');

% design_Method control
if design_Method ~= 1 && design_Method ~= 2 && design_Method ~= 3
    disp('Please check design method')
    return
end

starting_time=cputime;
%% 3. Call Optimization Algorithm
% Read data
Data
Input
clear dp_SI_M dp_US_M ffc_SI ffc_US ffci_SI ffci_US IX nbr_clms_sec
clear nbr_clms_tend nbr_rows_sec nbr_rows_tend nofdesignv nu_M
clear sec1 sec2 sec3 sec4 sec5 sec6 sec7 sec8 sec9 sec10 sec11 sec12...
    sec13 sec14 sec15 sec16 sec17 sec18 sec_all Section_Name Section_Name_Sorted...
    sort_sec_SI_M sort_sec_US_M
clear Tendon_SI_M_sorted Tendon_US_M Tendon_US_M_sorted

%Number of run
global ub lb n

N_run = 10;
GlobalBest = 1e40;
GlobalParam = zeros(n,1);
write = 0;
for i = 1:N_run
    % Optimization
	[Best, Param, hisx, hisy]= Optimization(ub,lb,i);
    
    fprintf('Run#: %i CurrentBest = %7.3e GlobalBest =%7.3e\n',i,Best,GlobalBest)
    if Best < GlobalBest
        GlobalBest = Best;
        GlobalParam = Param;
        Hisx = hisx;
        Hisy = hisy;
    end
end

%% Displaying results
write = write + 1;

global nu fp_c gama_conc Variable_section_properties Constant_section_properties Deflected_section_properties
 if GlobalBest > 1e39
    disp('Cannot find optimum design')
    [GlobalBest, const, fitness] = Computation(Param);
    return
 else
    disp('---------------------------------------------------------');
    disp('<strong>Optimum Results</strong>');
    fprintf('              Type or grade of concrete = %i\n', fp_c)
    fprintf('Ratio of instantaneous prestress losses = %7.4f\n', nu)
    fprintf('            Self weight of the concrete = %7.4f (kN/m^3)\n',gama_conc)
    [GlobalBest, const, fitness] = Computation(GlobalParam);
    fprintf('             Total Penalty = <strong>%7.4f\n</strong>',const)
    fprintf('                   fitness = <strong>%7.4f\n</strong>',fitness)
    disp(' ')
    fprintf('              Optimum Cost = <strong>%7.4f($)\n</strong>',GlobalBest)    
 end
 
disp('---------------------------------------------------------');
if design_Method == 1
    disp(Variable_section_properties)
elseif design_Method == 2
    disp(Constant_section_properties)
else
    disp(Deflected_section_properties)
end
disp('---------------------------------------------------------');

disp('Search History')
n = size(Hisx,2);
for i = 1:n
    fprintf('iter: %i GlobalBest: %6.3f\n', Hisx(i),  Hisy(i))
end

disp('------------------------------------------')
elapsed_time = (cputime - starting_time);
disp(['Elapsed time = ',num2str(elapsed_time),' seconds'])
disp(mfilename)
disp('Completed running')
disp(datestr(now))
disp('------------------------------------------')

%% Graphic
figure(3)
 plot(Hisx,Hisy,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)

xlabel('Iteration')
ylabel('Cost($)')
title('Search History of the Optimization Algorithm')

%% Saving Search History
if Opt_Method == 1
    if design_Method == 1
        % Delete existing file
        if isfile('ABC Variable Search History.png')
            delete('ABC Variable Search History.png')
        end
        % Create new file
        saveas(figure(3),'ABC Variable Search History.png')
        AVx = Hisx; AVy = Hisy;
        save ABC_Variable_history.mat AVx AVy
    elseif design_Method == 2
        % Delete existing file
        if isfile('ABC Constant Search History.png')
            delete('ABC Constant Search History.png')
        end
        % Create new file
        saveas(figure(3),'ABC Constant Search History.png')
        ACx = Hisx; ACy = Hisy;
        save ABC_Constant_history.mat ACx ACy
    else
        % Delete existing file
        if isfile('ABC Deflected Search History.png')
            delete('ABC Deflected Search History.png')
        end
        % Create new file
        saveas(figure(3),'ABC Deflected Search History.png')
        ADx = Hisx; ADy = Hisy;
        save ABC_Deflected_history.mat ADx ADy
    end
    
elseif Opt_Method == 2
    if design_Method == 1
        % Delete existing file
        if isfile('BSO Variable Search History.png')
            delete('BSO Variable Search History.png')
        end
        % Create new file
        saveas(figure(3),'BSO Variable Search History.png')
        BVx = Hisx; BVy = Hisy;
        save BSO_Variable_history.mat BVx BVy
    elseif design_Method == 2
        % Delete existing file
        if isfile('BSO Constant Search History.png')
            delete('BSO Constant Search History.png')
        end
        % Create new file
        saveas(figure(3),'BSO Constant Search History.png')
        BCx = Hisx; BCy = Hisy;
        save BSO_Constant_history.mat BCx BCy
    else
        % Delete existing file
        if isfile('BSO Deflected Search History.png')
            delete('BSO Deflected Search History.png')
        end
        % Create new file
        saveas(figure(3),'BSO Deflected Search History.png')
        BDx = Hisx; BDy = Hisy;
        save BSO_Deflected_history.mat BDx BDy
    end
    
else
    if design_Method == 1
        % Delete existing file
        if isfile('SSO Variable Search History.png')
            delete('SSO Variable Search History.png')
        end
        % Create new file
        saveas(figure(3),'SSO Variable Search History.png')
        SVx = Hisx; SVy = Hisy;
        save SSO_Variable_history.mat SVx SVy
    elseif design_Method == 2
         % Delete existing file
        if isfile('SSO Constant Search History.png')
            delete('SSO Constant Search History.png')
        end
        % Create new file
        saveas(figure(3),'SSO Constant Search History.png')
        SCx = Hisx; SCy = Hisy;
        save SSO_Constant_history.mat SCx SCy
    else
        % Delete existing file
        if isfile('SSO Deflected Search History.png')
            delete('SSO Deflected Search History.png')
        end
        % Create new file
        saveas(figure(3),'SSO Deflected Search History.png')
        SDx = Hisx; SDy = Hisy;
        save SSO_Deflected_history.mat SDx SDy
    end
end
