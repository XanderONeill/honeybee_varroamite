close all
clearvars

N_G = 5; N_T = 5; 

N_R_vector = 1:250;
drift_vector = 0.035:0.005:0.13;       %min=0.035, av=0.083, mx=0.13

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
j = 1;
for i = 1:length(drift_vector)
    drift = drift_vector(i)

    while Lower_limit(1, i) == 0 

        N_R = N_R_vector(j);         SP_p = [drift, 1000];   

        if N_R == 100
            break
        end

        N = N_G + N_T + N_R; % number of hives

        types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
        tspan = 0:1:10*30;

        options = odeset('NonNegative', 1:3*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        y0 = repmat([6, 14225, 10], 1, N);
        [~,y] = ode45(@(t,y) HoneyBee_space(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

        if y(end, 3*(N_G+N_T+N_R)-1) >= 10 
            Lower_limit(1, i) = N_R;
        else
            j = j+1
        end

        y = [];
    end
end

%% Mid Limit
for i = 1:length(drift_vector)
    drift = drift_vector(i)
    j = Lower_limit(i);

    while Mid_limit(1, i) == 0
               
        N_R = N_R_vector(j);         SP_p = [drift, 1000];     
        if N_R == 250
            break
        end
      
        N = N_G + N_T + N_R; % number of hives
        
        types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
        tspan = 0:1:10*20;
        
        options = odeset('NonNegative', 1:3*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        y0 = repmat([6, 14225, 10], 1, N);
        [~,y] = ode45(@(t,y) HoneyBee_space(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

        if y(end, 3*(N_G+N_T)-1) >= 10 
            Mid_limit(1, i) = N_R;
        else
            j = j+1
        end

        y = [];
    end
end

%% Upper Limit
% j = 50;
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
        tspan = 0:1:10*20;
        
        options = odeset('NonNegative', 1:3*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        y0 = repmat([6, 14225, 10], 1, N);
        [~,y] = ode45(@(t,y) HoneyBee_space(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

        if y(end, 3*(N_G)-1) >= 10 
            Upper_limit(1, i) = N_R;
        else
            j = j+1
        end

        y = [];
    end
end

save('Gen_tol_res_limits')
%% Figures 
close all 
clearvars
load('Gen_tol_res_limits.mat')

x = drift_vector;
y1 = Lower_limit./(N_T + N_G + Lower_limit);
y2 = Mid_limit./(N_T + N_G + Mid_limit);
y3 = Upper_limit./(N_T + N_G + Upper_limit);

fig = figure;
%subplot(2,1,1)
plot(x, y1,  'k-', 'LineWidth', 2)
hold on
plot(x, y2,  'k:', 'LineWidth', 2)
plot(x, y3,  'k-.', 'LineWidth', 2)

% patch([x fliplr(x)], [zeros(1,length(x)) fliplr(y1)], [0.6350 0.0780 0.1840], 'FaceAlpha', .2)
% patch([x fliplr(x)], [y1 fliplr(y2)], [0 0.4470 0.7410], 'FaceAlpha', .2)
% % patch([x fliplr(x)], [y1 fliplr(y3)], [0.9290 0.6940 0.1250], 'FaceAlpha', .2)
% patch([x fliplr(x)], [y2 ones(1,length(x))*100], [0.4660 0.6740 0.1880], 'FaceAlpha', .2)

xlim([0.035 0.13]);     xticks([0.04 0.07 0.1 0.13]);
xlabel('Drift parameter, \nu')

ylim([0 1]);            yticks([0 0.25 0.5 0.75 1]);
ylabel({'Proportion of resistant to total hives'});

set(gca, 'position', [0.15 0.12 0.8 0.65]);   
ax = gca;       ax.FontSize = 11;       set(gca,'box','off');
grid on;                    ax.GridColor = [0.5 .5 .5]; 
grid minor;
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;

