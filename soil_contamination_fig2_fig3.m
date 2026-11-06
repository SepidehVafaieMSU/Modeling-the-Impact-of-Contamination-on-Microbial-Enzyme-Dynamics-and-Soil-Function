%% soil_contamination_all_figures.m
% Reproduces every panel of Figure 2 and Figure 3 in
% Vafaie, Goodey, Krumins & Forgoston, "Modeling the Impact of Contamination
% on Microbial-Enzyme Dynamics and Soil Function" (bioRxiv 2026).
%
% Runs 4 scenarios once:  {uncontaminated, contaminated} x {vegetated, barren}
% and makes one figure window per panel (same as in the paper), saved to ./figures.
%
% State vector (same order in every run):
%   y = [C; U; D; M; E; CO2]
%   C = CSOC, U = USOC, D = DOC, M = microbial biomass, E = enzymes, CO2 = cumulative CO2
%
% The uncontaminated model is the SAME system with I_C = 0 and C(0) = 0:
% then C stays exactly 0, so r_ep = r_ep0 (Eq. 5), which is exactly the
% uncontaminated model of Section 2.1 (MICtoSOC in the old code = f1).
%
% The simulation is stopped automatically (event function) when C/U reaches 1,
% i.e. when r_ep = r_ep0*(1 - C/U) hits zero (Section 2.1).

clear; clc; close all;

%% ---------------- Parameters (Table 1) ----------------
P.T      = 20;          % temperature (deg C)
P.R      = 0.008314;    % kJ mol^-1 K^-1
P.rd     = 0.0002;      % microbial death rate (h^-1)
P.rep0   = 0.000005;    % enzyme production rate, no contamination (h^-1)
P.rel    = 0.001;       % enzyme loss rate (h^-1)
P.f1     = 0.5;         % fraction of dead biomass -> USOC
P.f2     = 0;           % fraction of dead biomass -> CSOC

% Assimilation (DOC uptake), Eqs. 7-8
P.kcatA = 1e8 * exp(-47 / (P.R*(P.T+273)));
P.KmA   = 0.1 + 0.01*P.T;
% USOC decomposition, Eqs. 16-17
P.kcatU = 1e8 * exp(-47 / (P.R*(P.T+273)));
P.KmU   = 500 + 5*P.T;
% CSOC decomposition, Eqs. 21-22 (k_catC = k_catU / 100)
P.kcatC = (1e8/100) * exp(-47 / (P.R*(P.T+273)));
P.KmC   = 500 + 5*P.T;
% Carbon use efficiency, Eq. 9
P.CUE   = 0.63 - 0.016*P.T;

tEnd_h = 30*8760;                       % 30 years in hours
y0_con = [0; 100; 0.5; 0.5; 0.01; 0];   % [C U D M E CO2] (C(0)=0 for both models)

Iveg = 0.0005;           % I_U = I_D in vegetated soil
Ibar = 0.0005/100;       % I_U = I_D in barren soil (100x lower)
IC   = 0.0005;           % press contamination input I_C

%% ---------------- Scenarios ----------------
%            name                    I_C   I_U    I_D
S = struct('name', {'Uncontaminated, vegetated', 'Uncontaminated, barren', ...
                    'Contaminated, vegetated',   'Contaminated, barren'}, ...
           'IC',   {0,    0,    IC,   IC}, ...
           'IU',   {Iveg, Ibar, Iveg, Ibar}, ...
           'ID',   {Iveg, Ibar, Iveg, Ibar});

opts      = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'Events',@stopWhenCequalsU);
optsShort = odeset('RelTol',1e-8, 'AbsTol',1e-10);

for k = 1:numel(S)
    rhs = @(t,y) soilModel(t, y, P, S(k).IC, S(k).IU, S(k).ID);

    % Long run: 30 years, or until C/U = 1
    [T, Y, te] = ode45(rhs, [0 tEnd_h], y0_con, opts);
    S(k).t_yr = T/8760;
    S(k).Y    = Y;
    S(k).stopped = ~isempty(te);

    % Short run for the DOC inset: first 10 hours, finely resolved
    tH = linspace(0, 10, 500);
    [~, Yh] = ode45(rhs, tH, y0_con, optsShort);
    S(k).t_h = tH;
    S(k).D_h = Yh(:,3);
end

%% ---------------- Print final values (compare with Table 2) ----------------
fprintf('\n%-28s %8s %9s %9s %10s %8s %10s %9s\n', ...
        'Scenario','t_end(yr)','C','U','D','M','E','CO2');
for k = 1:numel(S)
    y = S(k).Y(end,:);
    fprintf('%-28s %8.4f %9.3f %9.3f %10.7f %8.3f %10.3e %9.3f\n', ...
            S(k).name, S(k).t_yr(end), y);
end
fprintf('\n');

%% ---------------- Plot style ----------------
blue   = [0      0.4470 0.7410];   % vegetated
orange = [0.8500 0.3250 0.0980];   % barren
LW = 4;  axFS = 25;  labFS = 16;


UV = 1; UB = 2; CV = 3; CB = 4;     % scenario indices
yl = {'C (mg\cdotcm^{-3})', 'U (mg\cdotcm^{-3})', 'D (mg\cdotcm^{-3})', ...
      'M (mg\cdotcm^{-3})', 'E (mg\cdotcm^{-3})', 'CO_2 (mg\cdotcm^{-3})'};

% Panel list:  file name, state column, veg scenario, barren scenario, ylim, inset?
panels = {
  'Fig2a_CSOC_contaminated',  1, CV, CB, [0 150],   false
  'Fig2b_USOC_uncontaminated',2, UV, UB, [80 200],  false
  'Fig2c_USOC_contaminated',  2, CV, CB, [80 200],  false
  'Fig2d_DOC_uncontaminated', 3, UV, UB, [0 0.5],   true
  'Fig2e_DOC_contaminated',   3, CV, CB, [0 0.5],   true
  'Fig3a_M_uncontaminated',   4, UV, UB, [0 2.5],   false
  'Fig3b_M_contaminated',     4, CV, CB, [0 2.5],   false
  'Fig3c_E_uncontaminated',   5, UV, UB, [0 0.012], false
  'Fig3d_E_contaminated',     5, CV, CB, [0 0.012], false
  'Fig3e_CO2_uncontaminated', 6, UV, UB, [0 300],   false
  'Fig3f_CO2_contaminated',   6, CV, CB, [0 300],   false
};

%% ---------------- Draw every panel ----------------
for p = 1:size(panels,1)
    [fname, col, kv, kb, ylims, hasInset] = panels{p,:};

    fig = figure('Name', fname, 'Color', 'w', 'Position', [100 100 700 560]);
    ax  = axes(fig); hold(ax, 'on'); box(ax, 'on');

    plot(ax, S(kv).t_yr, S(kv).Y(:,col), '-', 'Color', blue,   'LineWidth', LW);
    plot(ax, S(kb).t_yr, S(kb).Y(:,col), ':', 'Color', orange, 'LineWidth', LW);

    xlim(ax, [0 30]); xticks(ax, 0:5:30); ylim(ax, ylims);
    set(ax, 'FontSize', axFS);
    xlabel(ax, 'Time (years)', 'FontSize', labFS);
    ylabel(ax, yl{col},       'FontSize', labFS);

    if hasInset   % zoom on DOC during the first 10 hours
        axIn = axes(fig, 'Position', [0.4 0.4 0.35 0.35]); hold(axIn,'on'); box(axIn,'on');
        plot(axIn, S(kv).t_h, S(kv).D_h, '-', 'Color', blue,   'LineWidth', LW);
        plot(axIn, S(kb).t_h, S(kb).D_h, ':', 'Color', orange, 'LineWidth', LW);
        xlim(axIn, [0 10]); ylim(axIn, [0 0.5]);
        set(axIn, 'FontSize', 15);
        xlabel(axIn, 'Time (hours)');
        ylabel(axIn, 'D (mg\cdotcm^{-3})');
    end

   
end


%% ======================= Local functions =======================
function dy = soilModel(~, y, P, IC, IU, ID)
    C = y(1); U = y(2); D = y(3); M = y(4); E = y(5);

    assim  = P.kcatA * M * D / (P.KmA + D);   % DOC assimilation, Eq. 2
    death  = P.rd * M;                        % Eq. 3
    rep    = P.rep0 * (1 - C/U);              % Eq. 5
    eprod  = rep * M;                         % Eq. 4
    eloss  = P.rel * E;                       % Eq. 12
    decU   = P.kcatU * E * U / (P.KmU + U);   % Eq. 15
    decC   = P.kcatC * E * C / (P.KmC + C);   % Eq. 20

    dC   = IC + P.f2*death - decC;                                        % Eq. 23
    dU   = IU + P.f1*death - decU;                                        % Eq. 18
    dD   = ID + (1-P.f1-P.f2)*death + decU + decC + eloss - assim;        % Eq. 25
    dM   = assim*P.CUE - death - eprod;                                   % Eq. 6
    dE   = eprod - eloss;                                                 % Eq. 13
    dCO2 = assim*(1 - P.CUE);                                             % Eq. 10

    dy = [dC; dU; dD; dM; dE; dCO2];
end

function [value, isterminal, direction] = stopWhenCequalsU(~, y)
    % r_ep = r_ep0*(1 - C/U) reaches 0 when C = U
    value      = y(2) - y(1);   % U - C
    isterminal = 1;             % stop the integration
    direction  = -1;            % only when U - C goes from + to -
end