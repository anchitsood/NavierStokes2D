
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
u_top = 1;

% Bottom wall of channel
stream_func_bottom = 5; % some user defined value
u_bottom = 0;

% Left wall of channel
stream_func_left = 5; % some user defined value
u_left = 0;

% Right wall of channel
stream_func_right = 5; % some user defined value
u_right = 0;

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

stream_func = zeros(elements_X+1,elements_Y+1);

stream_func(1, :) = stream_func_top;  
stream_func(:, 1) = stream_func_left;
stream_func(:, (elements_X + 1)) = stream_func_right;
stream_func((elements_Y + 1), :) = stream_func_bottom; 



%%% Step 2
%%% Start with vorticity defined everywhere (how? assume zero at all points for now) at t = 0. define now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vorticity = zeros(elements_X+1,elements_Y+1);


%%% Step 3
%%% use these 2 infos to define streamfunc everywhere at t = 0. define now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:elements_X + 1 - 1 
    for j = 2:elements_Y + 1 - 1
        sum1 = (delta_Y.^2 * (stream_func(i+1, j) + stream_func(i-1, j)));
        sum2 = (delta_X.^2 * (stream_func(i,j+1) + stream_func(i,j-1)));
        sum3 = (vorticity(i,j)*delta_X.^2 * delta_Y.^2);
        sum4 = sum1 + sum2 + sum3;
        stream_func(i,j) =  sum4  / (2 * (delta_Y.^2 + delta_X.^2));       
        %stream_func(i,j) = ((delta_Y.^2 * (stream_func(i+1, j) + stream_func(i-1, j))) + (delta_X.^2 * (stream_func(i,j+1) + stream_func(i,j-1))) + (vorticity(i,j)*delta_X.^2 * delta_Y.^2)) / (2 * (delta_Y.^2 + delta_X.^2));       
    end
end



%%% Step 4
%%% Increment time, recalculate vorticity values at the wall. taylor series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Boundary Condition for Top Wall
vorticity(:,1) = ((stream_func(:,1) - stream_func(:,2))* 2/(delta_Y).^2) + 2*u_top/(delta_Y);

%Boundary Condition for Bottom Wall
vorticity(:,elements_Y + 1) = ((stream_func(:,elements_Y + 1) - stream_func(:,elements_y + 1 - 1))* 2/(delta_Y).^2) + 2*u_bottom/(delta_Y);

%Boundary Condition for Left Wall
vorticity(1,:) = ((stream_func(1,:) - stream_func(2,:))* 2/(delta_X).^2) + 2*u_left/(delta_X);

%Boundary Condition for Right Wall
vorticity(elements_X + 1,:) = ((stream_func(elements_X + 1,:) - stream_func(elements_X+1-1,:))* 2/(delta_X).^2) + 2*u_right/(delta_X);





%%% Step 5
%%% recalculate vorticity values everywhere. navier stokes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clipped_vorticity = vorticity(2:m-1,2:n-1);
clipped_stream_func = stream_func(2:m-1,2:n-1);

vort_row_diffs = zeros(m - 2,n);
% vort_row_sums = zeros(m - 2,n);
vort_col_diffs = zeros(m,n - 2);
% vort_col_sums = zeros(m,n - 2);

strm_row_diffs = zeros(m - 2,n);
strm_row_sums = zeros(m - 2,n);
strm_col_diffs = zeros(m,n - 2);
strm_col_sums = zeros(m,n - 2);

for i = 1:(m - 2)
    vort_row_diffs(i,:) = vorticity(i + 2,:) - vorticity(i,:);
    % vort_row_sums(i,:) = vorticity(i + 2,:) + vorticity(i,:);
    strm_row_diffs(i,:) = stream_func(i + 2,:) - stream_func(i,:);
    strm_row_sums(i,:) = stream_func(i + 2,:) + stream_func(i,:);
end

for j = 1:(n - 2)
    vort_col_diffs = vorticity(:,j + 2) - vorticity(:,j);
    % vort_col_sums = vorticity(:,j + 2) + vorticity(:,j);
    strm_col_diffs(i,:) = stream_func(:,j + 2) - stream_func(:,j);
    strm_col_sums(i,:) = stream_func(:,j + 2) + stream_func(:,j);
end

vort_row_diffs = vort_row_diffs(:,2:n-1);
vort_col_diffs = vort_col_diffs(2:m-1,:);
% vort_row_sums = vort_row_sums(:,2:n-1);
% vort_col_sums = vort_col_sums(2:m-1,:);

strm_row_diffs = vort_row_diffs(:,2:n-1);
strm_col_diffs = vort_col_diffs(2:m-1,:);
strm_row_sums = vort_row_sums(:,2:n-1);
strm_col_sums = vort_col_sums(2:m-1,:);



term1 = - ((strm_row_diffs .* vort_col_diffs)/(4*delta_X*delta_Y));
term2 = ((strm_col_diffs .* vort_row_diffs)/(4*delta_X*delta_Y));
term3 = ((((delta_Y^2)*strm_row_sums) + ((delta_X^2)*strm_col_sums) - (2*clipped_stream_func*((delta_X^2) + (delta_Y^2))))/(Reynolds_no*(delta_X^2)*(delta_Y^2)));

grand_diff = term1 + term2 + term3;
vorticity_temp = (grand_diff * delta_T) + clipped_vorticity;

vorticity(2:m-1,2:n-1) = vorticity_temp;

% clear clipped_vorticity clipped_stream_func vort_row_diffs vort_col_diffs strm_row_diffs strm_row_sums strm_col_diffs strm_col_sums;
% clear term1 term2 term3 grand_diff vorticity_temp;



%%% Step 6
%%% Recalculate stream function everywhere. because boundary stream function doesn't change in time. continuity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 2:elements_X + 1 - 1 
    for j = 2:elements_Y + 1 - 1
        sum1 = (delta_Y.^2 * (stream_func(i+1, j) + stream_func(i-1, j)));
        sum2 = (delta_X.^2 * (stream_func(i,j+1) + stream_func(i,j-1)));
        sum3 = (vorticity(i,j)*delta_X.^2 * delta_Y.^2);
        sum4 = sum1 + sum2 + sum3;
        stream_func(i,j) =  sum4  / (2 * (delta_Y.^2 + delta_X.^2));       
        %stream_func(i,j) = ((delta_Y.^2 * (stream_func(i+1, j) + stream_func(i-1, j))) + (delta_X.^2 * (stream_func(i,j+1) + stream_func(i,j-1))) + (vorticity(i,j)*delta_X.^2 * delta_Y.^2)) / (2 * (delta_Y.^2 + delta_X.^2));       
    end
end










%%% Increment time, rinse and repeat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%% Step 7
%%% Calculate corresponding velocities in X and Y direction given stream function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros((elements_X + 1), (elements_Y + 1));
v = zeros((elements_X + 1), (elements_Y + 1));

for i =2: elements_X
    for j = 2:elements_Y
        u(i,j) = (stream_func(i,j+1)-stream_func(i,j-1))/2*delta_y;
        v(i,j) = (stream_func(i-1,j)-stream_func(i+1,j))/2*delta_x;
    end
end

%%
%%% Step 8
%%% Calculate corresponding Pressure in X and Y direction given stream function
pressure = zeros((elements_X + 1), (elements_Y + 1));


