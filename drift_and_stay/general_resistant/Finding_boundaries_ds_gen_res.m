close all
clearvars

N_G = 0; N_T = 20; 

N_R_vector = 1:250;
drift_vector = (0.035:0.005:0.13)/0.37;        %min=0.035, av=0.083, mx=0.13

Lower_limit = zeros(1, length(drift_vector));
Mid_limit = zeros(1, length(drift_vector));
Upper_limit = zeros(1, length(drift_vector));

%% Parameters uninfluenced by tolerant, resistant, regular honeybee scenario
        
%%% ----------- Brood Population Parameters ----------- %%%

a_H = 65500;        % Maximum birth rate coefficient
omega = 15200;      % Weigthing term for birth under Holling type-III
b_B = 0.97;         % Natural death rate of brood
m = 30/21;          % brood progression to HB

%%% ----------- Honeybee Population Parameters ----------- %%%
q_H = 1/70000;      % density-dependent death term for hosts
b_H = 30/14;

%%% ----------- Varroa Mite Parameters ----------- %%%
b_V = 0.81;         % natural death rate of VM
alpha_B = 0.2;     % maximum death rate of brood due to varroa
sigma = 8050;       % Half saturation on brood death from varroa

a_V_G = 0.000105;   % rate of new VM per average brood bee (general and tolerant)
a_V_R = 0.000094;   % rate of new VM per average brood bee (resistant)   
alpha_H_G = 3;      % mortality of honeybees due to varroa burden (general and resistant)
alpha_H_T = 1;      % mortality of honeybees due to varroa burden (tolerant)

%%% New Baseline
% baseline with a_H(t) 3/9,1,3/9,0
% Normal:       a_V = 0.0000105,     alpha_H = 3
% Resistant:    a_V = 0.0000094,     alpha_H = 3
% Tolerant:     a_V = 0.0000105,     alpha_H = 1 

BB_p = [a_H, omega, b_B, m];      % Bee Brood Params
HB_p = [q_H, b_H];                % Honeybee params
VM_p = [b_V, alpha_B, sigma];     % Varroa Params

type_params = [[a_V_G, alpha_H_G]; [a_V_R, alpha_H_T]];

%% Lower Limit
% how many resistant hives are required for the persistence of a 
% resistant population

for i = 1:length(drift_vector)
% for i = 1:3
    drift = drift_vector(i)

    for j = 1:30 

        N_R = N_R_vector(j);         
        SP_p = [drift, 1000];   

        N = N_G + N_T + N_R; % number of hives

        types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
        tspan = 0:1:12*20;

        options = odeset('NonNegative', 1:5*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        y0 = repmat([6, 14225, 10, 0, 0], 1, N);
        [~,y] = ode45(@(t,y) HoneyBee_disperse_and_stay(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

        if (y(end, 5*(N_G+N_T+N_R)-3) <= 10) || (isnan(y(end, 5*(N_G+N_T+N_R)-3)))
            Lower_limit(1, i) = N_R+1
        else
            break
        end

        y = [];
    end
end

% save('Lower_lims_tol_res')
%% Mid Limit - for Tolerant and Resistant Only
% % how many resistant hives are required for the persistence of both
% % tolerant and resistant hives
% 
% for i = 1:length(drift_vector)
% % for i = 3:length(drift_vector)
%     drift = drift_vector(i)
%     j = Lower_limit(i);
% 
%     while Mid_limit(1, i) == 0
% 
%         N_R = N_R_vector(j);         SP_p = [drift, 1000];     
% 
%         N = N_G + N_T + N_R; % number of hives
% 
%         types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
%         tspan = 0:1:10*40;
% 
%         options = odeset('NonNegative', 1:5*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
%         y0 = repmat([6, 14225, 10, 0, 0], 1, N);
%         [~,y] = ode45(@(t,y) HoneyBee_disperse_and_stay(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);
% 
%         if y(end, 5*(N_G+N_T)-3) >= 10 
%             Mid_limit(1, i) = N_R
%         else
%             j = j+1
%         end
% 
%         y = [];
%     end
% end
% 
% save('Mid_limit_Tol_res');


%% Upper Limit - for general and resistant combination
% how many resistant hives are required for the persistence of both
% general and resistant hives

j = 1;
for i = length(drift_vector):-1:1
% for i = 3:length(drift_vector)
    drift = drift_vector(i)

    while Upper_limit(1, i) == 0

        N_R = N_R_vector(j);         SP_p = [drift, 1000];     
        if N_R == 250
            break
        end

        N = N_G + N_T + N_R; % number of hives

        types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
        tspan = 0:1:10*40;

        options = odeset('NonNegative', 1:5*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        y0 = repmat([6, 14225, 10, 0, 0], 1, N);
        [~,y] = ode45(@(t,y) HoneyBee_disperse_and_stay(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

        if y(end, 5*(N_G)-3) >= 10 
            Upper_limit(1, i) = N_R
        else
            j = j+1
        end

        y = [];
    end
end

save('Gen_res')
%% Figures 
close all 
clearvars
load('Gen_res')

x = drift_vector;
y1 = Lower_limit./(N_T + N_G + Lower_limit);
% y2 = Mid_limit./(N_T + N_G + Mid_limit);
y2 = Upper_limit./(N_T + N_G + Upper_limit);

fig = figure;
plot(x, y1,  'k-', 'LineWidth', 2)
hold on
plot(x, y2,  'k-.', 'LineWidth', 2)

xlim([0.1 0.35]);     xticks([0.1 0.15 0.2 0.25 0.3 0.35]);
xlabel('Drifting rate, \eta')

ylim([0 1]);            yticks([0 0.25 0.5 0.75 1]);
ylabel({'Proportion of resistant hives'});

set(gca, 'position', [0.15 0.12 0.8 0.65]);   
ax = gca;       ax.FontSize = 11;       set(gca,'box','off');
grid on;                    ax.GridColor = [0.5 .5 .5]; 
grid minor;
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;

