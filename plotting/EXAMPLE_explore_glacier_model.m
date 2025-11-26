% OPTIONS - change as desired

% Decide which model to load
root = '/Users/badgeley/projects/seasonality/Cryo23/';

fric = 3;
indep = 5;
catchment = 36;
regularization = 100;
velocity_cutoff = 200;
rheology_option = 0;
start_time = 2018.5;
end_time = 2020.1;

model_name = ['Model_NW_TransientInversionRun_indep' num2str(indep) ...
              '_' num2str(start_time) '_' num2str(end_time) ...
              '_fric' num2str(fric) '_catchment' num2str(catchment) ...
              '_reg' num2str(regularization) ...
              'vel' num2str(velocity_cutoff) ...
              'prior0rhe' num2str(rheology_option)]; 

% Do you wish to load the model with only seasonal ice fronts?
load_frontal_model = 1;

% Specify times of interest (using decimal years)
t_start = 2019;
t_end = 2020;
t_peak = 2019.5;
t_trough = 2019.7;


%% SETUP FOR MODEL DOMAIN
if catchment == 35
    xlims = [-3.56e5,-3.27e5];
    ylims = [-1.59e6,-1.56e6];
elseif catchment == 36
    xlims = [-3.6e5,-3.3e5];
    ylims = [-1.53e6,-1.5e6];
elseif catchment == 37
    xlims = [-3.75e5,-3.4e5];
    ylims = [-1.51e6,-1.49e6];
elseif catchment == 48
    xlims = [-5e5,-4e5];
    ylims = [-1.275e6,-1.22e6];
elseif catchment == 23
    xlims = [-2.65e5,-2e5];
    ylims = [-1.9455e6,-1.89e6];
end

%% LOAD MODELS and VARIABLES

% load main model
md = loadmodel([root model_name]);

% Extract observations from model
[t_obs,v_obs] = extract_obs_from_model(md);

% Extract constants from model
g = md.constants.g;
rho_ice = md.materials.rho_ice;
rho_water = md.materials.rho_water;
yts = md.constants.yts;
X = md.mesh.x;
Y = md.mesh.y;

% Extract model results
%time (decimal year)
t = [md.results.TransientSolution(:).time];
%time indices
t_start_pos = find(abs(t-t_start)==min(abs(t-t_start)));
t_end_pos = find(abs(t-t_end)==min(abs(t-t_end)));
t_peak_pos = find(abs(t-t_peak)==min(abs(t-t_peak)));
t_trough_pos = find(abs(t-t_trough)==min(abs(t-t_trough)));

%mask
mask = [md.results.TransientSolution(:).MaskIceLevelset];
%make binary (1 is ice or land covered in a thin layer of ice, 0 is ocean)
mask(mask>0)=0;
mask(mask<0)=1;

%velocity (m/yr) and components
v = [md.results.TransientSolution(:).Vel];
vx = [md.results.TransientSolution(:).Vx];
vy = [md.results.TransientSolution(:).Vy];

%thickness (m)
H = [md.results.TransientSolution(:).Thickness];

%surface elevation (m)
S = [md.results.TransientSolution(:).Surface];

%base of ice (m) (S-H)
b = [md.results.TransientSolution(:).Base];

%bed elevation (m)
B = md.geometry.bed;

%basal friction coefficient
if fric==1 || fric==3
    C = [md.results.TransientSolution(:).FrictionCoefficient];
elseif fric==2
    C = [md.results.TransientSolution(:).FrictionC];
end

%basal shear stress alpha multiplier
alpha = [md.results.TransientSolution(:).FrictionAlpha2];

%basal shear stress (Pa) and components
taub = basal_shear_stress_alpha(alpha,v/yts);
taubx = basal_shear_stress_alpha(alpha,vx/yts);
tauby = basal_shear_stress_alpha(alpha,vy/yts);

%effective pressure (Pa) using model assumption of full connectivity to
%ocean
if fric==1 || fric==3
    N = effective_pressure_full_connectivity(H,b,g,rho_ice,rho_water,0);
elseif fric==2
    N = effective_pressure_full_connectivity(H,b,g,rho_ice,rho_water,0.01);
end

%surface slope (degrees) and components
[sx,sy,s] = surface_slope(md,S,2);

%driving stress (Pa) and components
taud = driving_stress(H,g,rho_ice,s);
taudx = driving_stress(H,g,rho_ice,sx);
taudy = driving_stress(H,g,rho_ice,sy);

%ice overburden pressure (Pa)
P_ice = ice_pressure(H,g,rho_ice);

%estimate CN from Budd sliding law
m = 6; %this can be anything - it doesn't matter the original friction law
CN = CN_from_Budd(taub,v./yts,m,m);

%estimate static C assuming water pressure equals zero at one time during
%the year
C_est = C_estimate_no_water(CN(:,t_start_pos:t_end_pos),...
                            P_ice(:,t_start_pos:t_end_pos));

%estimate evolving effective pressure
N_est = CN./C_est;

%estimate evolving water pressure
%you may get a warning of negative water pressure because C_est was
%estimated using only the part of the time series that is of interest
P_water_est = water_pressure_from_CN(CN,P_ice,C_est);

% load a model with only the seasonal frontal forcing if requested
% it will have the same friction coefficeint as your model of interest
% specified with an "_f" at the end of the variable
% NOTE: not all variables are loaded, but you can change that
if load_frontal_model
    modelFront_name = ['Model_NW_TransientInversionRun_indep0_2007_2022_fric' num2str(fric) '_catchment' num2str(catchment) '.mat'];
    mdf = loadmodel([root modelFront_name]);
    t_f = [mdf.results.TransientSolution(:).time];
    v_f = [mdf.results.TransientSolution(:).Vel];
end

%% INITIAL SPATIAL PLOT
plotmodel(md,'data',v(:,t_end_pos),...
             'mask',mask(:,t_end_pos), ...
             'xlim#all',xlims, ...
             'ylim#all',ylims,...
             'title','Velocity (m/yr)',...
             'figure',1);

%% PLOTS

% choose point location for time series plots
figure(1);
pos = getXYIndicesFromSpatialPlot(X,Y);

% redo spatial velocity figure with location plotted
plotmodel(md,'data',v(:,t_end_pos),...
             'mask',mask(:,t_end_pos), ...
             'xlim#all',xlims, ...
             'ylim#all',ylims,...
             'title','Velocity (m/yr)',...
             'figure',1);
hold on;
plot(X(pos),Y(pos),'ko','MarkerSize',10,'linewidth',2);
hold off;

% plot velocity time series for 'pos' location
figure(2);
plot(t,v(pos,:),'linewidth',2);
hold on;
plot(t_obs,v_obs(pos,:),'ko','linewidth',2);
if load_frontal_model
    plot(t_f,v_f(pos,:),'--','linewidth',2);
    legend({'Full seasonal', 'Obs', 'Frontal seasonal'},'Location','northwest');
else
    legend({'Frontal seasonal', 'Obs'},'Location','northwest');
end
xlim([t_start,t_end]);
xlabel('year');
ylabel('velocity (m/yr)');
title('Velocity')
hold off;

% plot taub spatial
plotmodel(md,'data',taub(:,t_end_pos)*(1e-6),'mask',mask(:,t_end_pos), ...
             'xlim#all',xlims, ...
             'ylim#all',ylims,...
             'caxis#all',[0,.2],...
             'title','Basal Shear Stress (MPa)',...
             'figure',3);
hold on;
plot(X(pos),Y(pos),'ko','MarkerSize',10,'linewidth',2);
hold off;

% plot components of taub and taud for 'pos' location
figure(4);
plot(t,taub(pos,:)*(1e-6),'k-','linewidth',2);
hold on;
plot(t,taud(pos,:)*(1e-6),'k--','linewidth',2);
plot(t,taubx(pos,:)*(1e-6),'b-','linewidth',2);
plot(t,taudx(pos,:)*(1e-6),'b--','linewidth',2);
plot(t,tauby(pos,:)*(1e-6),'r-','linewidth',2);
plot(t,taudy(pos,:)*(1e-6),'r--','linewidth',2);
legend({'tau_b','tau_d','tau_b^x','tau_d^x','tau_b^y','tau_d^y'},...
       'location','southwest');
xlabel('year');
ylabel('stress (MPa)');
xlim([t_start,t_end]);
title('Basal Shear and Driving Stresses');
hold off;

% Compare taub and taud at beginning and middle of year
plotmodel(md,'data',taub(:,t_start_pos)-taud(:,t_start_pos),...
             'data',taub(:,t_peak_pos)-taud(:,t_peak_pos),...
             'data',taub(:,t_trough_pos)-taud(:,t_trough_pos),...
             'data',taub(:,t_end_pos)-taud(:,t_end_pos),...
             'mask#1',mask(:,t_start_pos),...
             'mask#2',mask(:,t_peak_pos),...
             'mask#3',mask(:,t_trough_pos),...
             'mask#4',mask(:,t_end_pos),...
             'caxis#all',[-1e4,1e4],...
             'ylim#all',ylims,'xlim#all',xlims,...
             'title#1',['tau_b - tau_d: ' num2str(t_start)],...
             'title#2',['tau_b - tau_d: ' num2str(t_peak)],...
             'title#3',['tau_b - tau_d: ' num2str(t_trough)],...
             'title#4',['tau_b - tau_d: ' num2str(t_end)],...
             'figure',5);
hold on;
plot(X(pos),Y(pos),'ko','MarkerSize',10,'linewidth',2);
axs = findobj(gcf,'type','axes');
for ii = 2:length(axs)
    copyobj(axs(1).Children(1),axs(ii))
end
hold off;

% plot water and ice pressures for 'pos' location
figure(6); 
plot(t,P_water_est(pos,:),'linewidth',2); 
hold on; 
plot(t,P_ice(pos,:),'linewidth',2); 
xlabel('year');
ylabel('pressure (Pa)');
xlim([t_start,t_end]);
legend({'water pressure','ice pressure'});
title('Water (estimated) and Ice Pressures');
hold off;

% Water pressure change over time
plotmodel(md,'data',P_water_est(:,t_peak_pos)-P_water_est(:,t_start_pos),...
             'data',P_water_est(:,t_trough_pos)-P_water_est(:,t_peak_pos),...
             'data',P_water_est(:,t_end_pos)-P_water_est(:,t_trough_pos),...
             'mask#1',mask(:,t_peak_pos),...
             'mask#2',mask(:,t_trough_pos),...
             'mask#3',mask(:,t_end_pos),...
             'ylim#all',ylims,'xlim#all',xlims,...
             'caxis#all',[-5e6,5e6],...
             'title','Estimated water pressure change (Pa)',...
             'figure',7);
hold on;
plot(X(pos),Y(pos),'ko','MarkerSize',10,'linewidth',2);
axs = findobj(gcf,'type','axes');
for ii = 2:length(axs)
    copyobj(axs(1).Children(1),axs(ii))
end
hold off;