function [Sh, Ih, Rh] = simulate_absir(n_population,Iv0, T, infection_rate, recovery_rate,pod_size,gathering_chance)
% Simulate agent-based transmission model. Uses a graph to represent social
% connectivity between agents.
% 
% Note: You may find it helpful to summarize the state history matrices via
% summation. For instance sum(Ih, 1) will return the total number of
% infected persons at each timestep in the simulation.
%
% Inputs
%   n_population: Population size
%   Iv0 (column vector): Initial infection state
%   T (integer): Number of timesteps to simulate
%   infection_rate (float): Infection rate (probabilistic)
%   recovery_rate (float): Recovery rate (probabilistic)
% 
% Returns
%   Sh (matrix): Susceptible state history
%   Ih (matrix): Infected state history
%   Rh (matrix): Recovered state history

    % Setup
      % Dimensions of initial state
    dim = n_population;
    Ih = zeros(dim, T); % Infection history
    Ih(:, 1) = Iv0; % Record the initial state to infection history
    Rh = zeros(dim, T); % Recovery history
    
    % Construct an "action" helper function
    function [I, R,n_population,pod_size,gathering_chance] = action(I, R,n_population,pod_size,gathering_chance)
        %Remove random set of people
        remove = rand(1,n_population)>gathering_chance;
        IR = [I,R];
        store = IR(remove,:);
        IR(remove,:) = [];
        I = IR(:,1);
        R = IR(:,2);
        M = pod_maker(size(IR,1),pod_size);
        graph_M = graph(M);
        %figure();
        %plot(graph_M,'NodeColor','k', 'LineWidth',0.1)

        % Compute infection probabilities based on the social graph
        v_eff = infection_rate * M * I;
        v_pr = v_eff ./ (1 + v_eff);
        
        % Draw random values
        v_infect = rand(size(IR,1), 1) <= v_pr;
        v_recover = rand(size(IR,1), 1) <= recovery_rate;
        
        % Infect non-recovered individuals
        I = I | v_infect & (~R);
        % Recover infected individuals
        R = R | (I & v_recover);
        I = I & (~R);

        %recover people outside of simulation
        v_recover2 = rand(size(store,1), 1) <= recovery_rate;
        store(:,2) = store(:,2) | (store(:,1) & v_recover2);
        
        %merge people back together
        for g=1:n_population
            if remove(g) ==1
                newR(g) = store(1,2);
                newI(g) = store(1,1);
                store(1) = [];
            else
                newR(g) = R(1);
                newI(g) = I(1);
                R(1) = [];
                I(1) = [];
            end
        end
        I = newI;
        R = newR;
        n_population = n_population;
        pod_size = pod_size;
        gathering_chance = gathering_chance;
    end
        
    % Run simulation
    for i = 2:T
        [Ih(:, i), Rh(:, i),n_population,pod_size,gathering_chance] = action(Ih(:, i-1), Rh(:, i-1),n_population,pod_size,gathering_chance);
    end
    
    % Compute susceptible history
    Sh = ones(dim, T) - Ih - Rh;
end

