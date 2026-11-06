%% soil_contamination_fig6_fig7.m
% Reproduces Figure 6 and Figure 7 of
% Vafaie, Goodey, Krumins & Forgoston, "Modeling the Impact of Contamination
% on Microbial-Enzyme Dynamics and Soil Function" (bioRxiv 2026).
%
% Section 3.3: discrete (pulse) vs continuous (press) contamination,
% in vegetated and barren environments. One run; two windows (Fig. 6 and Fig. 7),
% each 3x2 like the paper:
%   left column  = discrete (pulse),   right column = continuous (press)
%   blue solid   = vegetated,          orange dotted = barren
%
%   Press: C(t=0) = 0 and I_C = 0.0001 mg cm^-3 h^-1 for 30 years.
%   Pulse: the same total amount (0.0001 * 30*8760 = 26.28 mg cm^-3) is
%          deposited at t = 0, i.e. C(t=0) = 26.28 and I_C = 0 afterwards.
%
% All variable names follow Table 1 of the paper:
%
%   State variables (mg cm^-3)            y = [C; U; D; M; E; CO2]
%     C    contaminated soil organic carbon (CSOC)
%     U    uncontaminated soil organic carbon (USOC)
%     D    dissolved organic carbon (DOC)
%     M    microbial biomass
%     E    extracellular enzymes
%     CO2  cumulative CO2 respired
%
%   Time t is in hours (30 years = 30*8760 h); plots show years.
%
% A run is stopped when C = U, i.e. when r_ep = r_ep0*(1 - C/U) reaches 0
% (this does not happen in these scenarios; all runs reach 30 years).

clear; clc; close all;

%% ================= Table 1: parameters =================
P.R      = 0.008314;   % ideal gas constant (kJ mol^-1 K^-1)
P.T      = 20;         % temperature (deg C)

P.r_d    = 0.0002;     % microbial death rate (h^-1)
P.r_ep0  = 0.000005;   % enzyme production rate without contamination (h^-1)
P.r_el   = 0.001;      % enzyme loss rate (h^-1)
P.f1     = 0.5;        % fraction of dead microbial biomass added to USOC
P.f2     = 0;          % fraction of dead microbial biomass added to CSOC

% Carbon use efficiency (Eq. 9)
P.CUE_0  = 0.63;       % CUE at 0 deg C
P.CUE_s  = -0.016;     % slope (deg C^-1)

% DOC assimilation (Eqs. 7-8)
P.kcat_A0 = 1e8;       % pre-exponential coefficient (h^-1)
P.Ea_A    = 47;        % activation energy (kJ mol^-1)
P.Km_A0   = 0.1;       % Km_A at 0 deg C (mg cm^-3)
P.Km_As   = 0.01;      % slope (mg cm^-3 degC^-1)

% USOC decomposition (Eqs. 16-17)
P.kcat_U0 = 1e8;
P.Ea_U    = 47;
P.Km_U0   = 500;
P.Km_Us   = 5;

% CSOC decomposition (Eqs. 21-22); k_catC = k_catU / 100
P.kcat_C0 = 1e6;
P.Ea_C    = 47;
P.Km_C0   = 500;
P.Km_Cs   = 5;

% ----- calculated parameters (Table 1, "Calculated") -----
P.CUE    = P.CUE_0 + P.CUE_s*P.T;                          % Eq. 9  -> 0.31
P.kcat_A = P.kcat_A0 * exp(-P.Ea_A / (P.R*(P.T+273)));     % Eq. 7  -> 0.4176
P.Km_A   = P.Km_A0 + P.Km_As*P.T;                          % Eq. 8  -> 0.3
P.kcat_U = P.kcat_U0 * exp(-P.Ea_U / (P.R*(P.T+273)));     % Eq. 16 -> 0.4176
P.Km_U   = P.Km_U0 + P.Km_Us*P.T;                          % Eq. 17 -> 600
P.kcat_C = P.kcat_C0 * exp(-P.Ea_C / (P.R*(P.T+273)));     % Eq. 21 -> 0.004176
P.Km_C   = P.Km_C0 + P.Km_Cs*P.T;                          % Eq. 22 -> 600

%% ================= Table 1: initial values and time =================
% C_0 depends on the contamination mode (set below)
U_0   = 100;    % USOC
D_0   = 0.5;    % DOC
M_0   = 0.5;    % microbial biomass
E_0   = 0.01;   % extracellular enzymes
CO2_0 = 0;      % CO2

t_end = 30*8760;            % 30 years (h)

%% ================= External inputs and contamination modes =================
I_env   = [0.0005, 0.0005/100];   % I_U = I_D:  1 = vegetated, 2 = barren (100x lower)
envName = {'Vegetated', 'Barren'};

I_C_press = 0.0001;               % press contamination rate (mg cm^-3 h^-1)
C_pulse   = I_C_press * t_end;    % = 26.28 mg cm^-3, same total amount as press

modeName = {'Discrete (pulse)', 'Continuous (press)'};
%            C_0       I_C
modes = { C_pulse,  0;            % 1 = pulse
          0,        I_C_press };  % 2 = press

%% ================= Run all 4 simulations =================
opts      = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'Events',@stopWhenCequalsU);
optsShort = odeset('RelTol',1e-8, 'AbsTol',1e-10);
t_hours   = linspace(0, 10, 500);   % first 10 hours, for the DOC insets

R = struct();
for m = 1:2                               % contamination mode
    C_0 = modes{m,1};
    I_C = modes{m,2};
    y_0 = [C_0; U_0; D_0; M_0; E_0; CO2_0];
    for e = 1:2                           % environment
        I_U = I_env(e);  I_D = I_env(e);
        rhs = @(t,y) soilModel(t, y, P, I_C, I_U, I_D);

        [t, y] = ode45(rhs, [0 t_end], y_0, opts);
        R(m,e).t_yr = t/8760;
        R(m,e).y    = y;

        [~, yh] = ode45(rhs, t_hours, y_0, optsShort);
        R(m,e).D_hours = yh(:,3);
    end
end

%% ================= Final values (compare with Table 4) =================
fprintf('\n%-20s %-10s %9s %9s %9s %10s %7s %11s %9s\n', ...
        'Mode','Env','t_end(yr)','C','U','D','M','E','CO2');
for m = 1:2
    for e = 1:2
        fprintf('%-20s %-10s %9.4f %9.3f %9.3f %10.7f %7.3f %11.4e %9.3f\n', ...
                modeName{m}, envName{e}, R(m,e).t_yr(end), R(m,e).y(end,:));
    end
end
fprintf('\n');

%% ================= Plots =================
blue   = [0      0.4470 0.7410];   % vegetated
orange = [0.8500 0.3250 0.0980];   % barren
cols   = {blue, orange};
style  = {'-', ':'};               % vegetated solid, barren dotted
LW = 2.5;  FS = 13;
ylab = {'C (mg\cdotcm^{-3})', 'U (mg\cdotcm^{-3})', 'D (mg\cdotcm^{-3})', ...
        'M (mg\cdotcm^{-3})', 'E (mg\cdotcm^{-3})', 'CO_2 (mg\cdotcm^{-3})'};
letters = 'abcdef';

% rows of each figure: {state column, ylim}; columns: pulse / press
fig6rows = {1, [0 30];  2, [90 140]; 3, [0 0.5]};    % C, U, D
fig7rows = {4, [0 2.5]; 5, [0 0.01]; 6, [0 250]};    % M, E, CO2
figs      = {fig6rows, fig7rows};
figTitles = {'Figure 6', 'Figure 7'};

for f = 1:2
    rows = figs{f};
    fig = figure('Name', figTitles{f}, 'Color', 'w', 'Position', [60+40*f 40 1000 1100]);
    for r = 1:3
        col = rows{r,1};  ylims = rows{r,2};
        for m = 1:2                                   % pulse / press
            ax = subplot(3, 2, 2*(r-1)+m); hold(ax,'on'); box(ax,'on');
            for e = 1:2                               % vegetated / barren
                plot(ax, R(m,e).t_yr, R(m,e).y(:,col), style{e}, ...
                     'Color', cols{e}, 'LineWidth', LW);
            end
            xlim(ax, [0 30]); xticks(ax, 0:5:30); ylim(ax, ylims);
            set(ax, 'FontSize', FS);
            xlabel(ax, 'Time (years)'); ylabel(ax, ylab{col});
            title(ax, sprintf('(%s) %s', letters(2*(r-1)+m), modeName{m}), ...
                  'FontWeight', 'normal');

            % ----- insets -----
            if col == 1 && m == 1          % Fig 6(a): C under pulse, zoomed
                axIn = makeInset(fig, ax, [0.35 0.15 0.55 0.45]);
                for e = 1:2
                    plot(axIn, R(m,e).t_yr, R(m,e).y(:,1), style{e}, ...
                         'Color', cols{e}, 'LineWidth', 1.5);
                end
                xlim(axIn, [0 30]); ylim(axIn, [25.9 26.3]);
                xlabel(axIn, 'Time (years)'); ylabel(axIn, ylab{1});
            elseif col == 3                % Fig 6(e,f): D, first 10 hours
                axIn = makeInset(fig, ax, [0.35 0.35 0.55 0.55]);
                for e = 1:2
                    plot(axIn, t_hours, R(m,e).D_hours, style{e}, ...
                         'Color', cols{e}, 'LineWidth', 1.5);
                end
                xlim(axIn, [0 10]); ylim(axIn, [0 0.5]);
                xlabel(axIn, 'Time (hours)'); ylabel(axIn, ylab{3});
            end
        end
    end
    lg = legend(subplot(3,2,2), envName, 'Location', 'northwest', 'FontSize', 10);
    lg.Box = 'off';
end

%% ================= Local functions =================
function axIn = makeInset(fig, ax, rel)
    % rel = [left bottom width height] as fractions of the parent axes
    p = ax.Position;
    axIn = axes(fig, 'Position', [p(1)+rel(1)*p(3), p(2)+rel(2)*p(4), rel(3)*p(3), rel(4)*p(4)]);
    hold(axIn, 'on'); box(axIn, 'on');
    set(axIn, 'FontSize', 8);
end

%% ================= Model (Eqs. 2-25) =================
function dydt = soilModel(~, y, P, I_C, I_U, I_D)
    C = y(1); U = y(2); D = y(3); M = y(4); E = y(5);

    DOC_assimilation  = P.kcat_A * M * D / (P.Km_A + D);   % Eq. 2
    microbial_death   = P.r_d * M;                         % Eq. 3
    r_ep              = P.r_ep0 - P.r_ep0 * (C/U);         % Eq. 5
    enzyme_production = r_ep * M;                          % Eq. 4
    enzyme_loss       = P.r_el * E;                        % Eq. 12
    USOC_decomp       = P.kcat_U * E * U / (P.Km_U + U);   % Eq. 15
    CSOC_decomp       = P.kcat_C * E * C / (P.Km_C + C);   % Eq. 20

    dC   = I_C + P.f2*microbial_death - CSOC_decomp;                  % Eq. 23
    dU   = I_U + P.f1*microbial_death - USOC_decomp;                  % Eq. 18
    dD   = I_D + (1 - P.f1 - P.f2)*microbial_death + USOC_decomp ...
           + CSOC_decomp + enzyme_loss - DOC_assimilation;            % Eq. 25
    dM   = DOC_assimilation*P.CUE - microbial_death - enzyme_production;  % Eq. 6
    dE   = enzyme_production - enzyme_loss;                           % Eq. 13
    dCO2 = DOC_assimilation*(1 - P.CUE);                              % Eq. 10

    dydt = [dC; dU; dD; dM; dE; dCO2];
end

function [value, isterminal, direction] = stopWhenCequalsU(~, y)
    % r_ep = r_ep0*(1 - C/U) reaches 0 when C = U
    value      = y(2) - y(1);   % U - C
    isterminal = 1;
    direction  = -1;
end