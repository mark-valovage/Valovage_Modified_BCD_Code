function [ChangePoint_Indexes] = Get_ChangePoints_Windowed(mu0, kappa0, alpha0, beta0, lambda, hazard_func, X_entire, T, WINDOW_LENGTH)
%Get_ChangePoints_Windowed:  This function will take as input times series data and
%default starting values, and returns changepoint indexes (calculated 
% offline).  Code modeled after McKay's Change Detection code.  
%   Input variables (Refer to Adams'/MacKay's original paper).
%       Starting values
%       - mu0
%       - kappa0
%       - alpha0
%       - kappa0
%       Function characteristics
%       - lambda
%       - hazard_func
%       Time series values
%       - X_entire
%       - T
%       WINDOW_LENGTH (Size of the partition window to use).  
%
%   Output variables (Computed offline)
%       ChangePoint_Indexes
    
    %% Step 1:  Partition the time series arrays
    % Name change (for now, change in a bit...)
    window_start_index = 1;
    window_end_index = WINDOW_LENGTH;
    
    ChangePoints_Found = [];

    while(window_start_index <= size(X_entire,1))
        
        if (window_end_index > size(X_entire,1))
            window_end_index = size(X_entire,1);
        end
                         
        X = X_entire(window_start_index:window_end_index);
        T = size(X,1);
        
%         disp(['Window start index is ' num2str(window_start_index)])
%         disp(['Window end index is ' num2str(window_end_index)])
%         disp(['Running Change Detection on subwindow from ' num2str(window_start_index) ' to ' num2str(window_end_index)])

        %% Step 2:  Run Change Detection to Get Log Smears and Maxes
        %******Original Code from Adams/MacKay******
        R = zeros([T+1 T]);

        % At time t=1, we actually have complete knowledge about the run
        % length.  It is definitely zero.  See the paper for other possible
        % boundary conditions.
        R(1,1) = 1;

        %Initialize paramters
        muT    = mu0;
        kappaT = kappa0;
        alphaT = alpha0;
        betaT  = beta0;

        % Allocate space for maximums.
        maxes  = zeros([T+1]);

        % Loop over the data window like we're seeing it all for the first time.
        for t=1:T
          % Evaluate the predictive distribution for the new datum under each of
          % the parameters.  This is the standard thing from Bayesian inference.
          predprobs = studentpdf(X(t), muT, ...
                                 betaT.*(kappaT+1)./(alphaT.*kappaT), ...
                                 2 * alphaT);
  
          % Evaluate the hazard function for this interval.
          H = hazard_func([1:t]');
  
          % Evaluate the growth probabilities - shift the probabilities down and to
          % the right, scaled by the hazard function and the predictive
          % probabilities.
          R(2:t+1,t+1) = R(1:t,t) .* predprobs .* (1-H);
  
          % Evaluate the probability that there *was* a changepoint and we're
          % accumulating the mass back down at r = 0.
          R(1,t+1) = sum( R(1:t,t) .* predprobs .* H );
  
          % Renormalize the run length probabilities for improved numerical
          % stability.
          R(:,t+1) = R(:,t+1) ./ sum(R(:,t+1));

          % Update the parameter sets for each possible run length.
          muT0    = [ mu0    ; (kappaT.*muT + X(t)) ./ (kappaT+1) ];
          kappaT0 = [ kappa0 ; kappaT + 1 ];
          alphaT0 = [ alpha0 ; alphaT + 0.5 ];
          betaT0  = [ beta0  ; betaT + (kappaT .*(X(t)-muT).^2)./(2*(kappaT+1)) ];
          muT     = muT0;
          kappaT  = kappaT0;
          alphaT  = alphaT0;
          betaT   = betaT0;
  
          % Store the maximum.
          maxes(t) = find(R(:,t)==max(R(:,t)));
  
        end
        %******End Original Code from Adams/MacKay******
             
        %% Step 3:  Extract changepoints (offline)
        % Extract Changepoints from Maxes 
        diffs_maxes = diff(maxes(:,1));
        ChangeStates_End = [find(diffs_maxes ~= 1), diffs_maxes(find(diffs_maxes ~= 1))];

        offset = window_start_index - 1;
        End_Sequence_Index = ChangeStates_End(end,1);
        
        % Iteratively trace the active run sequence back to its origin.
        while (End_Sequence_Index > 0)
            Sequence_Length = maxes(End_Sequence_Index,1);
            Start_Sequence_Index = End_Sequence_Index - Sequence_Length;
            if (Start_Sequence_Index > 0)
                % This index is the start of an 'active' run sequence.
                % Save the time index as a changepoint.
                ChangePoints_Found = [Start_Sequence_Index + offset, Sequence_Length; ChangePoints_Found];
            end
            End_Sequence_Index = Start_Sequence_Index - 1;
        end

        % Move indexes to next window segment
        window_start_index = window_end_index + 1;
        window_end_index = window_end_index + WINDOW_LENGTH;
    end

    %% Step 4:  Return ChangePoints
    if(isempty(ChangePoints_Found))
        ChangePoint_Indexes = [];
    else
        ChangePoint_Indexes = sort(ChangePoints_Found(:,1));
    end

end

