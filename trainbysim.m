classdef trainbysim < handle

    properties
        %%%%%%%%%%%% ATTRIBUTES %%%%%%%%%%%%%%%
        reservoir_type
        Nr
        rho
        tau
        inputScaling
        biasScaling
        lambda
        connectivity
        Win
        Wb
        Wr
        randweight
        initialstate
        
        Wout
        Woutmat
        train_internalState
        train_internalStatedot
        train_reservoirReadout
        train_reservoirReadoutWashed
        train_reservoirTarget
        train_reservoirTargetWashed
        predict_internalState
        predict_internalStatedot
        predict_reservoirReadout
        sizeinput
        nodenuminput
        sizeoutput
        delaylen
        timestep
        ifxnorm
        xnorm
        delaynorm
        normmethod

    end

    methods
        function esn = trainbysim(Nr, varargin)
            %%%%%%%%%%%% CONSTRUCTOR CLASS %%%%%%%%%%%%%
            esn.Nr = Nr;  %number of nodes
            esn.tau = 10;  %timescale
            esn.inputScaling = 1;
            esn.biasScaling = 1;
            esn.lambda = 1;
            esn.connectivity = 1;
            esn.sizeinput = 1;
            esn.randweight = true;
            esn.nodenuminput = 1;
            esn.sizeoutput = 1;
            esn.initialstate = zeros(esn.Nr,1);
            esn.delaynorm = nan;
            esn.xnorm = 1;
            esn.train_reservoirReadoutWashed = [];
            esn.Woutmat = [];   %matrix for output weights

            numvarargs = length(varargin);
            for i = 1:2:numvarargs
                switch varargin{i}
                    case 'timeConst', esn.tau = varargin{i+1};
                    case 'spectralRadius', esn.rho = varargin{i+1};
                    case 'inputScaling', esn.inputScaling = varargin{i+1};
                    case 'biasScaling', esn.biasScaling = varargin{i+1};
                    case 'regularization', esn.lambda = varargin{i+1};
                    case 'connectivity', esn.connectivity = varargin{i+1};
                    case 'sizeinput', esn.sizeinput = varargin{i+1};
                    case 'nodenuminput', esn.nodenuminput = varargin{i+1};
                    case 'sizeoutput', esn.sizeoutput = varargin{i+1};
                    case 'randweight', esn.randweight = varargin{i+1};
                    case 'delaylen', esn.delaylen = varargin{i+1};
                    case 'timestep', esn.timestep = varargin{i+1};
                    case 'normmethod', esn.normmethod = varargin{i+1};
                    case 'ifxnorm', esn.ifxnorm = varargin{i+1};
                    case 'reservoir_type', esn.reservoir_type = varargin{i+1};

                    otherwise, error('the option does not exist');
                end
            end

            esn.train_reservoirTarget = cell(esn.sizeoutput,1);
            esn.train_reservoirTargetWashed = cell(esn.sizeoutput,1);
            for i = 1:esn.sizeoutput
                esn.train_reservoirTargetWashed{i} = [];
            end
            esn.Wout = cell(1,1);
            esn.Woutmat = zeros(1, esn.Nr+esn.sizeinput+1);  %Initial weights so training has no external forces

            switch esn.reservoir_type
                case 'ESN'
                    esn.Win = esn.inputScaling * ones(esn.Nr, 1);
                    esn.Win(2:2:end) = -esn.Win(2:2:end);            %alternatate at a everyother rate 
                    esn.Wb = esn.biasScaling * (rand(esn.Nr, 1) * 2 - 1);  
                    esn.Wr = full(sprand(esn.Nr, esn.Nr, esn.connectivity));   %generate random matrix with a size of Nr and a connectivity 
                    esn.Wr(esn.Wr ~= 0) = esn.Wr(esn.Wr ~= 0) * 2 - 1;
                    esn.Wr = esn.Wr * (esn.rho / max(abs(eig(esn.Wr))));   %scale the matrix so the spectral radius meet the desired.

                case {'FFNN'}
                    esn.Win = rand(esn.Nr, 1);
                    esn.Wb = esn.biasScaling * (rand(esn.Nr, 1) * 2 - 1);  
                    esn.Wr = zeros(esn.Nr, esn.Nr);
                    esn.Wr(2,1) = rand();
                    esn.Wr(3,2) = rand();
                    esn.Wr(4,3) = rand();



                    
        otherwise, error('This option does not exist')
            end
        end

 

        function [] = clearrecord(esn, varargin)
            %%%%%%RESET VALUES - MAINTAINS GOOD RESERVOIRS %%%%%%%%%%%
            numvarargs = length(varargin);
            for i = 1:2:numvarargs
                switch varargin{i}
                    case 'timeConst', esn.tau = varargin{i+1};
                    case 'inputScaling', esn.inputScaling = varargin{i+1};
                    case 'regularization', esn.lambda = varargin{i+1};
                    case 'sizeinput', esn.sizeinput = varargin{i+1};
                    case 'nodenuminput', esn.nodenuminput = varargin{i+1};
                    case 'sizeoutput', esn.sizeoutput = varargin{i+1};
                    case 'randweight', esn.randweight = varargin{i+1};
                    case 'delaylen', esn.delaylen = varargin{i+1};
                    case 'timestep', esn.timestep = varargin{i+1};
                    case 'normmethod', esn.normmethod = varargin{i+1};
                    case 'ifxnorm', esn.ifxnorm = varargin{i+1};
                                          
                    otherwise, error('the option does not exist');
                end
            end

            esn.train_internalState = [];
            esn.train_internalStatedot = [];
            esn.train_reservoirReadout = [];
            esn.train_reservoirReadoutWashed = [];
            esn.train_reservoirTarget = cell(esn.sizeoutput,1);
            esn.train_reservoirTargetWashed = cell(esn.sizeoutput,1);
            for i = 1:esn.sizeoutput
                esn.train_reservoirTargetWashed{i} = [];
            end
            esn.predict_internalState = [];
            esn.predict_internalStatedot = [];
            esn.predict_reservoirReadout = [];
            esn.xnorm = 0;
            esn.delaynorm = nan;
            esn.predict_internalState = [];
            esn.predict_reservoirReadout = [];
            esn.xnorm = ones(esn.sizeinput,1);
        end


        function [x, target] = traintest(esn, simname, stoptime)
            %%%%%%%%%%%%%%% NORMALISE INPUT VALUES %%%%%%%%%%%%%%%%%%%

            [const, x, target, internalState, internalStatedot] = esn.runsim(simname, stoptime);

            v1 = x(1:end/2,:);  %take the first half of values
            if esn.ifxnorm
                esn.xnorm = max(abs(v1),[],2)/2;    %adjust so the xnorm to force inputs to be unified
                esn.xnorm = [esn.xnorm; esn.xnorm];
            end
        end



         function [] = traindatacollect(esn, simname, washout, stoptime)
            %%%%%%%%%%%%%%%%%% COLLECT TRAINING DATA %%%%%%%%%%%%%%%%%%%%
            [const, x, target, internalState, internalStatedot] = esn.runsim(simname, stoptime);
            esn.train_internalState = internalState;
            esn.train_internalStatedot = internalStatedot;

            esn.train_reservoirReadout = [const'; x; esn.train_internalState'];   %Create vector of [constant; robot_states; reservoir_states] for future training 
            esn.train_reservoirReadoutWashed = [esn.train_reservoirReadoutWashed esn.train_reservoirReadout(:,washout+1:end)]; %Adjust so it does not take the beginning (wait for synchronization) values (washout is the amount of entries to skip) 
            
            for i = 1:esn.sizeoutput   %For future development - repeat this for all desired varaibles
                esn.train_reservoirTarget{i} = target{i};
                esn.train_reservoirTargetWashed{i} = [esn.train_reservoirTargetWashed{i}; esn.train_reservoirTarget{i}(washout+1:end)];
            end
         end

 
         function [controloutput] = robotcontrol(esn, simname, stoptime)       % Use weights to calculate reservoir readout torque
            %%%%%%%%%%%%%%%%%%%%%% CONTROL ROBOT %%%%%%%%%%%%%%%%%%%%%%%%%
            [const, x, target, internalState] = esn.runsim(simname, stoptime);

            controloutput = [const'; x; internalState'; target{1}'];    %Input out results from sim to [constant; robot_states; reservoir_states]
            
        end

        function [train_output, train_target] = train(esn)   %Calculate the reservoir weights based on the training data
            %%%%%%%%%%%%%%%%% CALCULATE TRAINING WEIGHTS %%%%%%%%%%%%%%%%%%%%%%%
            train_output = cell(esn.sizeoutput,1);
            train_target = cell(esn.sizeoutput,1);

            for i = 1:esn.sizeoutput    
                X = esn.train_reservoirReadoutWashed; 
                Y = esn.train_reservoirTargetWashed{i};
                disp("Size of X")
                disp(size(X))
                disp("Size of Y ")
                disp(size(Y))
                esn.Wout{i} = Y'*X'*inv(X*X'+esn.lambda*eye(size(X,1)));    %Linear regression formula
                esn.Woutmat = esn.Wout{i};
                disp("These are the weights")
                disp(esn.Woutmat)
                train_output{i} = esn.Wout{i}*esn.train_reservoirReadout;
                train_output{i} = train_output{i}';
                train_target{i} = esn.train_reservoirTarget{i};

            end
        end


        function [const, x, target, internalState, internalStatedot] = runsim(esn, simname, simTime)
            %%%%%%%%%%%%%%%%%%% RUN SIMULATION AND COLLECT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
            switch simname
                case "Pendulum_network_control"
                    out = sim(simname, 'StopTime', num2str(simTime));   %Call the Simulink model 

                    for i = 1:esn.Nr
                        internalState(:,i) = out.yout{2}.Values.S_esn.data(i,1,:); %Collect the states over time
                        internalStatedot(:,i) = out.yout{2}.Values.Sdot_esn.data(i,1,:);   %Collect the state derivatives over time 
                    end
                    
                    const = squeeze(out.yout{1}.Values.const.Data);  %collect bias over time
                    x = squeeze(out.yout{1}.Values.x.Data);  %collect robot states (theta, theta, theta*, thetadot*)
                    target = cell(1);
                    target{1} = squeeze(out.yout{1}.Values.target.data);  %acceleration (torque) over time

                otherwise 
                    fprintf('Simulation model: %s\n', simname);
                    error('the simulation does not exist');
            end
        end
    end
end