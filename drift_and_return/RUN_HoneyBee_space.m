close all
clearvars

% clearvars -except drift_vector Lower_limit Upper_limit N_T
N_G = 0; N_T = 20; N_R = 60;

drift = 0.083;       %min=0.035, av=0.083, mx=0.13

%% Parameter justification

% SEE ANDY'S ONENOTE DOCUMENT FOR MORE INFORMATION; will write this up into
% a pdf document.

% omega = 15200
%           Khoury et al. (2011) and Booton et al. (2017), but adjusted
%           into a Holling-type III functional form
% m = 30/21
%           lots of sources for this one - Khoury et al. (2013)

% K_H = 50000, K_B = 25000
%           maximum honeybee colony size. 
%           Indicated as this in Schmickl and Crailsheim (2007)

% a_H = 65500 
%           maximum honeybee reproduction per month. We initially had this 
%           as 3000 per day, but many papers seem to indicate a maximum of
%           approximately 2000 - citation:
%           Messan et al. (2021) indicates max of 2500 (simulation)
%           DeGrandiHoffman et al. ('89) indicates between 2000 and 3000 (simulation)
%           Sumpter and Martin ('04) indicates 1500 - from Snodgrass (1925) 
%           Schmeckl and Crailsheim ('07) indicates 2000 - Bodenheimer (1937)
%           Shpigler et al. (2022) cites Laidlaw and Page (1997) - 2000
%                   Paper title ^ - The influences of illumination on...      

% b_B = 0.97        
%           natural death rate of brood - 
% q_H = 1/70000     
%           density-dependent term for honeybee regulation

%           these terms are chosen such that, in the absence of varroa, the
%           brood and honeybee populations can reach the carrying 
%           capacities of K_B and K_H, respectively. 

% b_V = 0.81        
%           natural death rate of varroa - Messan et al. (2021)

% a_V = 0.000065   
%           taken to get a steady state density of Brood, in the presence 
%           of varroa, to be half the carrying capacity of Brood in the
%           absence of varroa. i.e. B* = K_B = 25000/2 = 12500

%           NEED TO FIT: to normal and resistant bee types
%           What is normal level of MPHB in resistant bees (assume ~1 for
%           now)

% alpha_B = 0.72
% sigma = 8050
%           alpha_B - maximum death rate of brood due to infestation of
%           varroa
%           sigma - half saturation for influence of brood density on brood
%           death from varroa
%           This term is obtained from Messan et al. (2021), along with
%           respective parameters.

% alpha_H   
%           mortality of honeybees due to varroa burden, 
% 
%           FITTED: to normal and tolerant bee types.
%           McMahon 2016 - normal last 3 years, tolerant last 4 years


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
SP_p = [drift, 1000];             % Dispersal params [nu, H_0];

% H_0 = 1000 seems reasonable - half saturation point
% nu = eta*p,         where:
%       p is the probability of a varroa moving from one bee to another
%       p = 0.37 from data (find Lewis source)

%% Model Simulations
N = N_G + N_T + N_R; % number of hives

types = [ones(N_G,1); 2*ones(N_T,1);    3*ones(N_R,1)];
type_params = [[a_V_G, alpha_H_G]; [a_V_R, alpha_H_T]];
tspan = 0:1:12*50;

options = odeset('NonNegative', 1:3*N, 'RelTol', 1e-9, 'AbsTol', 1e-9);
y0 = repmat([6, 14225, 10], 1, N);
[t,y] = ode45(@(t,y) HoneyBee_space(t, y, BB_p, HB_p, VM_p, SP_p, types, type_params), tspan, y0, options);

%% Plotting

close all

%% General Hives
if N_G > 0
fig1 = figure;

% ---- Honeybees and brood ----%
subplot(2,1,1)
plot(t(1:8*12+1), y(1:8*12+1,3*N_G-2), 'k-.', 'LineWidth', 1.5)
hold on
plot(t(1:8*12+1), y(1:8*12+1,3*N_G-1), '-', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250])

xlim([0 8*12]);     xticks([0 12 24 36 48 60 72 84 96]);       xticklabels('');
ylim([0 5.1e4]);    yticks([0 2e4 4e4]);    ylabel('Density per hive');

ax = gca;       ax.FontSize = 11;       set(gca,'box','off');
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
set(gca, 'position', [0.12 0.55 0.6 0.4])
leg1 = legend('Brood', 'Honeybee', 'Position', [0.82 0.8 0.1 0.1]);     
set(leg1,'Box','off');

% title(['T =' num2str(N_G), ', R=' num2str(N_R), ', drift=av'])

% ---- Varroa Plots ----%
subplot(2,1,2)
% ---- Varroa Densities ----%
yyaxis left
plot(t(1:8*12+1), y(1:8*12+1,3*N_G), '-', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])

xlim([0 8*12]);     xticks([0 12 24 36 48 60 72 84 96]);
ylim([0 3000]);     yticks([0 1e3 2e3 3e3]);    ylabel('Mites');
ax = gca;           ax.FontSize = 11;           ax.YColor = [0.6350 0.0780 0.1840];

% ---- Varroa per hundred bee ----%
yyaxis right
plot(t(1:8*12+1), 100*y(1:8*12+1,3*N_G)./y(1:8*12+1,3*N_G-1), '--','LineWidth', 1.5, 'Color', [0 0.4470 0.7410])

xlim([0 8*12]);     xticks([0 12 24 36 48 60 72 84 96]);     xlabel('Time, t (months)');
ylim([0 15]);       yticks([0 5 10 15]);                  ylabel('MPHB');

ax.FontSize = 11;           ax.YColor = [0 0.4470 0.7410];
set(gca,'box','off');       set(gca, 'position', [0.12 0.12 0.6 0.4]);       
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
leg2 = legend('Mites', 'MPHB', 'Position', [0.82 0.6 0.1 0.1]);  
set(leg2,'Box','off')
end
%% Tolerant Hives
if N_T > 0
fig2 = figure;

% ---- Honeybees and brood ----%
subplot(2,1,1)
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T)-2), 'k-.', 'LineWidth', 1.5)
hold on
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T)-1), '-', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);       xticklabels('');
ylim([0 5.1e4]);    yticks([0 2e4 4e4]);    ylabel('Density per hive');

ax = gca;       ax.FontSize = 11;       set(gca,'box','off');
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
set(gca, 'position', [0.12 0.55 0.6 0.4])
leg1 = legend('Brood', 'Honeybee', 'Position', [0.82 0.8 0.1 0.1]);     
set(leg1,'Box','off');

% title(['T =' num2str(N_G), ', R=' num2str(N_R), ', drift=av'])

% ---- Varroa Plots ----%
subplot(2,1,2)
% ---- Varroa Densities ----%
yyaxis left
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T)), '-', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);
ylim([0 9000]);     yticks([0 3e3 6e3 9e3]);    ylabel('Mites');
ax = gca;           ax.FontSize = 11;           ax.YColor = [0.6350 0.0780 0.1840];

% ---- Varroa per hundred bee ----%
yyaxis right
plot(t(1:10*12+1), 100*y(1:10*12+1,3*(N_G + N_T))./y(1:10*12+1,3*(N_G + N_T)-1), '--','LineWidth', 1.5, 'Color', [0 0.4470 0.7410])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);     xlabel('Time, t (months)');
xtickangle(0);
ylim([0 45]);       yticks([0 15 30  45]);                  ylabel('MPHB');

ax.FontSize = 11;           ax.YColor = [0 0.4470 0.7410];
set(gca,'box','off');       set(gca, 'position', [0.12 0.12 0.6 0.4]);       
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
leg2 = legend('Mites', 'MPHB', 'Position', [0.82 0.6 0.1 0.1]);  
set(leg2,'Box','off')
end

%% Resistant Hives
if N_R > 0
fig3 = figure;

% ---- Honeybees and brood ----%
subplot(2,1,1)
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T + N_R)-2), 'k-.', 'LineWidth', 1.5)
hold on
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T + N_R)-1), '-', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);       xticklabels('');
ylim([0 5.1e4]);    yticks([0 2e4 4e4]);    ylabel('Density per hive');

ax = gca;       ax.FontSize = 11;       set(gca,'box','off');
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
set(gca, 'position', [0.12 0.55 0.6 0.4])
leg1 = legend('Brood', 'Honeybee', 'Position', [0.82 0.8 0.1 0.1]);     
set(leg1,'Box','off');

% title(['T =' num2str(N_G), ', R=' num2str(N_R), ', drift=av'])

% ---- Varroa Plots ----%
subplot(2,1,2)
% ---- Varroa Densities ----%
yyaxis left
plot(t(1:10*12+1), y(1:10*12+1,3*(N_G + N_T + N_R)), '-', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
% hold on
% plot(t(1:10*12+1), y(end-(10*12):end,3*(N_G + N_T + N_R)), '-.','LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);
ylim([0 3000]);     yticks([0 1e3 2e3 3e3]);    ylabel('Mites');
ax = gca;           ax.FontSize = 11;           ax.YColor = [0.6350 0.0780 0.1840];

% ---- Varroa per hundred bee ----%
yyaxis right
plot(t(1:10*12+1), 100*y(1:10*12+1,3*(N_G + N_T + N_R))./y(1:10*12+1,3*(N_G + N_T + N_R)-1), '--','LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
% hold on
% plot(t(1:10*12+1), 100*y(end-(10*12):end,3*(N_G + N_T + N_R))./y(end-(10*12):end,3*(N_G + N_T + N_R)-1), ':','LineWidth', 1, 'Color', [0 0.4470 0.7410])

xlim([0 10*12]);     xticks([0 12 24 36 48 60 72 84 96 108 120]);     xlabel('Time, t (months)');
xtickangle(0);
ylim([0 9]);       yticks([0 3 6 9]);                  ylabel('MPHB');

ax.FontSize = 11;           ax.YColor = [0 0.4470 0.7410];
set(gca,'box','off');       set(gca, 'position', [0.12 0.12 0.6 0.4]);       
grid on;                    ax.GridColor = [0.5 .5 .5]; 
ax.GridLineStyle = '--';    ax.GridAlpha = 0.5;
leg2 = legend('Mites', 'MPHB', 'Position', [0.82 0.6 0.1 0.1]);  
% leg2 = legend('Mites', 'Mites (long time)','MPHB','MPHB (long time)', 'Position', [0.82 0.6 0.1 0.1]);  
set(leg2,'Box','off')
end