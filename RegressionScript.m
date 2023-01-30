%% Simulate a diffusion system
%  This system uses numerical differentation to simulate the diffusion of a
%  component out of a spherical bead into a container of fixed volume.
clf
clc
clear

%% Read data from file (commented, update if file changes)
N = 20; % ~, number of discrete nodes in finite difference

opts = detectImportOptions('DiffusionData.xlsx', 'Sheet', 'Diffusion data');
opts = setvartype(opts,{'HydrogelName'},'string');
RawData = readtable('DiffusionData.xlsx', opts);

Types = ["Glycerol", "Xylitol", "Sorbitol"];

for i = 1:length(Types)
    temp = RawData(RawData.HydrogelName == Types(i),:);
    
    % Parameters associated with this sphere
    Data.(Types(i)).p.V = temp.VolumeSolution_ml_(1) / 1e6; % mL -> m3, volume of container
    Data.(Types(i)).p.R = 0.5*temp.SphereDiameter_cm_(1) * 1e-2; % mm -> m, radius of sphere
    Data.(Types(i)).p.C0 = temp.InitialSphereConcentration_mM_(1); % mM -> mol/m3, initial concentration of component in bead
    Data.(Types(i)).p.A = 4*pi*Data.(Types(i)).p.R^2;    % m2, surface area of the sphere
    Data.(Types(i)).p.N = N;          % ~, number of discrete nodes in finite difference

    Data.(Types(i)).estPi = temp.PartitionCoefficient(1); % ~, partition coefficient estimate

    % Measured data for this sphere
    Data.(Types(i)).t = temp.Time_hours_ * 3600; % hrs -> s, measurement times
    Data.(Types(i)).cf = temp.SolutionConcentration_mM_; % mM -> mol/m3, solution concentration
end



%% Start simulation
clc
for i = 1:length(Types)
    Current = Data.(Types(i));

    
    D = 0.5e-9;       % m2/s, diffusion coefficient
    Pi = Current.estPi;  % ~, partition coefficient Pi = cf/c(r=R) = c_fluid / c inside bead surface
    
    
    opts = optimoptions('lsqnonlin','FiniteDifferenceType','central', ...
                        'TypicalX', [1e-11, 1], ...
                        'StepTolerance', 1e-12);
    
    p_regress = lsqnonlin(@(p_regress) CalculateSquaredDifference(p_regress, Current), [D, Pi], [0, 0],[], opts);
    Current.D = p_regress(1);
    Current.Pi = p_regress(2);
    
    
    %% Bootstrap
    figure(i)
    t = linspace(0, Current.t(end), 100);
    N = 5000;
    n = length(Current.t);
    for j = 1:N
        index = randsample(n, n, true);
        p_regress = lsqnonlin(@(p_regress) CalculateSquaredDifference(p_regress, Current, index), [D, Pi], [0, 0],[1e-7 10], opts);
        Current.Bootstrap.D(j) = p_regress(1);
        Current.Bootstrap.Pi(j) = p_regress(2);
        
        subplot(2,1,1)
        plot(t/3600, CalculateSolutionConcentration(Current.Bootstrap.D(j), Current.Bootstrap.Pi(j), Current.p, t), 'LineWidth',0.5, 'Color', 0.6*[1 1 1]);
        hold on
        sprintf('Type %s, Iteration %i / %i', Types(i), j, N)
    end
    
    subplot(2,1,1)
    plot(t/3600, CalculateSolutionConcentration(Current.D, Current.Pi, Current.p, t), 'b', ...
         Current.t/3600, Current.cf, 'ro', 'LineWidth', 2);
    title(sprintf('%s', Types(i)));
    xlabel('Time (hrs)'); ylabel('Solution concentration (mM)');

    subplot(2,2,3);
    histogram(Current.Bootstrap.D*1e9);
    title(sprintf('D = %0.3f +- %0.3f x 10^-^9 m^2/s', mean(Current.Bootstrap.D*1e9), std(Current.Bootstrap.D*1e9)))
    xlabel('Diffusion coefficient (m^2/s) x 10^9')
    subplot(2,2,4);
    histogram(Current.Bootstrap.Pi);
    title(sprintf('Pi = %0.3f +- %0.3f', mean(Current.Bootstrap.Pi), std(Current.Bootstrap.Pi)))
    xlabel('Partition coefficient')
    
    
    %% Identifiability analysis for D 
    figure(i+3)
    N = 50; % Number of segments in identifiability
    opts_Pi = optimoptions('lsqnonlin','FiniteDifferenceType','central', ...
                        'StepTolerance', 1e-12);
    DVec = logspace(-13, -4, N);
    for j = 1:length(DVec)
        [~, Current.D_SSE(j)] = lsqnonlin(@(Pi) CalculateSquaredDifference_Pi(Pi, DVec(j), Current), Current.Pi, [0], [], opts_Pi);
    end
    
    subplot(2,1,1);
    loglog(DVec, Current.D_SSE)
    xlabel('D'); ylabel('SSE');
    title(sprintf('Identifiability analysis for D (%s)', Types(i)))
    
    %% Identifiability analysis for Pi 
    opts_D = optimoptions('lsqnonlin','FiniteDifferenceType','central', ...
                        'TypicalX', 1e-11, ...
                        'StepTolerance', 1e-12);
    PiVec = logspace(-2, 2, N);
    for j = 1:length(PiVec)
        [~, Current.Pi_SSE(j)] = lsqnonlin(@(D) CalculateSquaredDifference_D(D, PiVec(j), Current), Current.D, [0], [], opts_D);
    end
    
    subplot(2,1,2);
    loglog(PiVec, Current.Pi_SSE)
    xlabel('\Pi'); ylabel('SSE');
    title(sprintf('Identifiability analysis for Pi (%s)', Types(i)))
    
    Data.(Types(i)) = Current;
end

save(['Data', datestr(datetime('now'), 'HHMMSS'),'.mat'], "Data");


%% Regression functions
function ErrorVector = CalculateSquaredDifference(p_regress, Current, index)
    D = p_regress(1);
    Pi = p_regress(2);

    cf = CalculateSolutionConcentration(D, Pi, Current.p, Current.t);

    ErrorVector = cf' - Current.cf;
    if nargin == 3
        ErrorVector = ErrorVector(index);
    end
end

function ErrorVector = CalculateSquaredDifference_Pi(Pi, D, Current)
    cf = CalculateSolutionConcentration(D, Pi, Current.p, Current.t);
    ErrorVector = cf' - Current.cf;
end

function ErrorVector = CalculateSquaredDifference_D(D, Pi, Current)
    cf = CalculateSolutionConcentration(D, Pi, Current.p, Current.t);
    ErrorVector = cf' - Current.cf;
end


%% Simulation functions
function [cf, c] = CalculateSolutionConcentration(D, Pi, p, t)

M = ConstructMatrix(Pi, p);

c0 = [p.C0*ones(p.N,1); 0];
c = zeros(p.N+1, length(t));
for i = 1:length(t)
    c(:,i) = expm(D*M*t(i))*c0;
end
cf = c(end,:);  % Concentration in solution

end

function M = ConstructMatrix(Pi, p)
    dr = p.R/(p.N+1);
    r = dr:dr:p.R;
    
    UpperDiagonal = diag(r(2:end)./r(1:end-1), 1);
    LowerDiagonal = [diag(r(1:end-2)./r(2:end-1), -1) zeros(p.N,1); zeros(1, p.N+1)];
    MainDiagonal = [-diag(2*ones(p.N,1)) zeros(p.N,1); zeros(1, p.N+1)];
    M = UpperDiagonal + LowerDiagonal + MainDiagonal;
    M(end-1,end) = M(end-1,end)/Pi;     % The last element in the vector will be cf = Pi*c(r=R)
    M(end, end-2:end) = -[((p.A/p.V)*(dr/2))  (-(p.A/p.V)*(4*dr/2))  ((p.A/p.V)*(3*dr/2)/Pi)];
    
    M = M/dr^2;
end