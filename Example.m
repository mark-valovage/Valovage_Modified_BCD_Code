clearvars;

% Add path to supporting functions
addpath('.\utils\');

% Load example data
load('BreadMaker.mat');

% Parameters for change detection. These are all domain agnostic and
% not tuned to any dataset.

% Hazard function 
%(Domain agnostic; see Adams/MacKay 2007 and original OBCD code)
lambda        = 200;
hazard_func   = @(r) constant_hazard(r, lambda);

% Starting statistical parameters 
% (Also domain agnostic; see Adams/MacKay 2007 and original OBCD code)
params.mu0    = 0;
params.kappa0 = 1;
params.alpha0 = 1;
params.beta0  = 1;

% Window length (see Section 3.3 of Valovage-AAMAS2017 paper)
params.WINDOW_LENGTH = 1000;

% Extract changepoint indexes using modified function
ChangePoint_Indexes = Get_ChangePoints_Windowed(params.mu0, ...
                                                params.kappa0, ...
                                                params.alpha0, ...
                                                params.beta0, ...
                                                lambda, ...
                                                hazard_func, ...
                                                BreadMaker.power, ...
                                                size(BreadMaker.power,1), ...
                                                params.WINDOW_LENGTH);


if (isempty(ChangePoint_Indexes))
    disp('No changePoints found')

else
    disp('ChangePoints found at time indexes:')
    disp(ChangePoint_Indexes)
end
                                            
% Plot time series
plot(BreadMaker.time, BreadMaker.power); grid; hold on;
xlim([0 201]);  % Specific for this example...

title('Breadmaker Example');

xlabel('Time (seconds)');
ylabel('Power (Watts)');


% Plot changepoints detected
ChangePoint_line_y_coords = [140,150];   % Specific for this example...
for i = 1:size(ChangePoint_Indexes,1)
    line([BreadMaker.time(ChangePoint_Indexes(i)),BreadMaker.time(ChangePoint_Indexes(i))],ChangePoint_line_y_coords,'Color','r', 'LineWidth',2);
end

legend('Breadmaker real power data stream', 'Changepoint');

hold off;