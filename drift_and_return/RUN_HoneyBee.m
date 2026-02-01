% close all
clearvars

RUN = 1; %'Regular';
% RUN = 2; %'Tolerant';
% RUN = 3; %'Resistant';

%%% ------------ Holling type-III birth function -------------%%%

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
%           approximately 2000 - we use 2000. citation:
%           Messan et al. (2021), from Sumpter and Martin 2007, indicates max of 1500 (simulation)
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

% a_V = 0.000061   
%           taken to get a steady state density of Brood, in the presence 
%           of varroa, to be half the carrying capacity of Brood in the
%           absence of varroa. i.e. B* = K_B = 25000/2 = 12500

%           NEED TO FIT: to normal and resistant bee types
%           What is normal level of MPHB in resistant bees (assume ~1 for
%           now)

% alpha_B = 0.72
%           Obtained from Messan et al. (2021) Figures
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
b_B = 0.97;         % Natural death rate of brood messan range
m = 30/21;          % brood progression to HB

%%% ----------- Honeybee Population Parameters ----------- %%%
q_H = 1/70000;      % density-dependent death term for hosts
b_H = 30/14;            % natural death of Honeybee

%%% ----------- Varroa Mite Parameters ----------- %%%
b_V = 0.81;         % natural death rate of VM
alpha_B = 0.2;   % maximum death rate of brood due to varroa
sigma = 8050;       % Half saturation on brood death from varroa

if RUN == 1
    a_V = 0.000105;    % rate of new VM per average brood bee
    alpha_H = 3;        % mortality of honeybees due to varroa burden
elseif RUN == 2
    a_V = 0.000105;    % rate of new VM per average brood bee
    alpha_H = 1;        % mortality of honeybees due to varroa burden
elseif RUN == 3
    a_V = 0.000094;    % rate of new VM per average brood bee
    alpha_H = 3;        % mortality of honeybees due to varroa burden
end

%%% New Baseline
% baseline with a_H(t) 3/9,1,3/9,0
% Normal:       a_V = 0.000105,     alpha_H = 3
% Resistant:    a_V = 0.000094,     alpha_H = 3
% Tolerant:     a_V = 0.000105,     alpha_H = 1 


% baseline with a_H(t) 8/9,1,6/9,1/9
% Normal:       a_V = 0.000061,     alpha_H = 3
% FITTING - DIE WITHIN APPROX 3 YEARS,  MONTHS 18: 0<MPHB<9 
%                               (WITH INTERQUARTILE RANGE OF 0<MPHB<2.5)
%                                       MONTHS 21: 0<MPHB<16 
%                               (WITH INTERQUARTILE RANGE OF 0<MBHPH<6.5)
% Resistant:    a_V < 0.000052,     alpha_H = 3.
% Tolerant:     a_V = 0.000061,     alpha_H = 1. 

BB_p = [a_H, omega, b_B, m];
HB_p = [q_H, b_H];
VM_p = [a_V, b_V, alpha_B, sigma, alpha_H];

%% Model Simulations

options = odeset('NonNegative', 1:3);%, 'RelTol', 1e-9, 'AbsTol', 1e-9);
y0 = [6, 14225, 10];    % initial conditions
tspan = 0:1:12*50;          

[t,y] = ode15s(@(t,y) HoneyBee(t, y, BB_p, HB_p, VM_p), tspan, y0, options);

%% Plotting
% close all

figure
% ---- Honeybees and brood ----%
subplot(2,1,1)
plot(t(1:5*12+1), y(1:5*12+1, 1),  'k-.', 'LineWidth', 1.5)
hold on
plot(t(1:5*12+1), y(1:5*12+1,2), '-', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250])
xlim([0 5*12])
xticks([0 12 24 36 48 60])
xticklabels('')
ylim([0 5.1e4])
yticks([0 2e4 4e4])
ylabel('Density per hive')
ax = gca;
ax.FontSize = 11;
% title('Regular')
set(gca,'box','off')
set(gca, 'position', [0.12 0.55 0.6 0.4])
leg1 = legend('Brood', 'Honeybee', 'Position', [0.82 0.8 0.1 0.1]);
set(leg1,'Box','off')
grid on
%grid minor
ax.GridColor = [0.5 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;


% ---- Varroa Plots ----%
subplot(2,1,2)

% ---- Varroa Densities ----%
yyaxis left
plot(t(1:5*12+1), y(1:5*12+1,3), '-', 'LineWidth', 1.5, 'Color', [0.6350 0.0780 0.1840])
hold on
if RUN == 3
    plot(t(1:5*12+1), y(end-(5*12):end,3), '-.','LineWidth', 1, 'Color', [0.6350 0.0780 0.1840])
end
xlim([0 5*12])
xticks([0 12 24 36 48 60])
if RUN == 1
    ylim([0 3000])
    yticks([0 1e3 2e3 3e3])
elseif RUN == 2
    ylim([0 8000])
    yticks([0 2e3 4e3 6e3])
elseif RUN == 3
    ylim([0 550])
    yticks([0 200 400])
end
ylabel('Mites')
ax = gca;
ax.FontSize = 11;
ax.YColor = [0.6350 0.0780 0.1840];

% ---- Varroa per hundred bee ----%
yyaxis right
plot(t(1:5*12+1), 100*y(1:5*12+1,3)./y(1:5*12+1,2), '--','LineWidth', 1.5, 'Color', [0 0.4470 0.7410])
hold on
if RUN == 3
    plot(t(1:5*12+1), 100*y(end-(5*12):end,3)./y(end-(5*12):end,2), ':','LineWidth', 1, 'Color', [0 0.4470 0.7410])
end
xlim([0 5*12])
xticks([0 12 24 36 48 60])
xlabel('Time, t (months)')
if RUN == 1
    ylim([0 15])
    yticks([0 5 10 15])
elseif RUN == 2
    ylim([0 40])
    yticks([0 10 20 30 40])
elseif RUN == 3
    ylim([0 2.75])
    yticks([0 1 2])
end
ylabel('MPHB')
if RUN == 3
    leg2 = legend('Mites', 'Mites (long time)', 'MPHB', 'MPHB (long time)', 'Position', [0.82 0.6 0.1 0.1]);
else
    leg2 = legend('Mites', 'MPHB', 'Position', [0.82 0.6 0.1 0.1]);
end
set(leg2,'Box','off')
ax = gca;
ax.FontSize = 11;
set(gca,'box','off')
set(gca, 'position', [0.12 0.12 0.6 0.4])
ax.YColor = [0 0.4470 0.7410];
grid on
%grid minor
ax.GridColor = [0.5 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

% figure (2)
% plot(t, y(:,1), 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250])
% hold on
% plot(t, y(:,2), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980])
% f = find(y(:,2)<10000);
% % xline(f(1)-1, 'LineWidth', 2)
% xlim([0 tspan(end)])
% xticks([0 36 72 108 144 180 216])
% xticklabels('')
% ylim([0 5.1e4])
% yticks([0 2.5e4 5e4])
% ylabel('Density per hive')
% ax = gca;
% ax.FontSize = 11;
% % title('Regular')
% set(gca,'box','off')
% set(gca, 'position', [0.15 0.7 0.5 0.25])
% leg1 = legend('Brood', 'Honeybee','DEAD', 'Position', [0.72 0.8 0.1 0.1]);
% set(leg1,'Box','off')
% grid on
% %grid minor
% ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;


