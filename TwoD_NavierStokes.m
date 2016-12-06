
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description

% Members: 
% Bla Bla
% Bla
% Bla


%% Begin code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem parameters

%%% Define element sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify problem domain by side lengths of the rectangular area
length_X = 10;
length_Y = 10;

% Specify number of elements to split into
elements_X = 10; % = n_X
elements_Y = 10; % = n_Y

% Specify element size
delta_X = length_X/elements_X;
delta_Y = length_Y/elements_Y;



%%% Define user specified quantities (boundary conditions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Top wall of channel
stream_func_top = 5; % some user defined value

% Bottom wall of channel
stream_func_bottom = 5; % some user defined value

% Left wall of channel
stream_func_left = 5; % some user defined value

% Right wall of channel
stream_func_right = 5; % some user defined value

% Stream function overall is stored in a matrix called stream_func
syms stream_func;



%%% Define time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify time step and number of time steps
delta_T = 1;
stepsof_T = 50; % number of time steps

% Or, specify total time and size of time step
total_T = 50;
delta_T = 1;

stepsof_T = total_T/delta_T;


%% Solution steps

%%%% Begin time loop

%%% Step 1
%%% Start with streamfunction defined on boundaries at t = 0. define now or
%%% in parameteres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Step 2
%%% Start with vorticity defined everywhere (how? assume zero at all points for now) at t = 0. define now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Step 3
%%% use these 2 infos to define streamfunc everywhere at t = 0. define now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








for timestep_count = 1:stepsof_T
%%% Step 4
%%% Increment time, recalculate vorticity values at the wall. taylor series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Step 5
%%% recalculate vorticity values everywhere. navier stokes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Step 6
%%% Recalculate stream function everywhere. because boundary stream function doesn't change in time. continuity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










end
%%% Increment time, rinse and repeat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




