%%%%%%%%%%%% ENVIRONMENT PARAMETERS %%%%%%%%%%%%
g = 9.80665;  %Gravity
robot_damping_coeff = 0; 
torquelim = 1; % 1 for fully acutated <0.08 for underactuated

%%%%%%%%%%% RESERVOIR PARAMETERS %%%%%%%%%%%%
Win = eye(2); %For initlisation of reservoir parameters 
SR = 1.1;  %Scale factor for largest 
timeConstant = 0.01;
normmethod = 1;
ifxnorm = 1;
samplingrate = 0.01;
stoptime = 20;
washout = 0.03*stoptime/samplingrate;  %Initial amount of indices to ignore in training data 
sizeoutput = 1;  %onl
sizeinput = 4;
amp_factor = 0.01;   %Scaling factor for xnorms


%%%%%%%%%% INITIALISATION VALUES %%%%%%%%%%%%%%%%
Kp = 10;   %proportional control coefficient 
Kd = 1;    %dervaitve contol coefficient
Ki = 10;   %intergal control coefficient 

desired_theta = pi;  %the desired angle of the pendulum when using PID controller

PID_on = false; %PID is a robot controller, setting PID_on as false means the pendulum drops freely, PID_on as true mean PID controller tries to control to desired angle
control_ref = [1; 0];  %Used for initial 
simTime = 10;



