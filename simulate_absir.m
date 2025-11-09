function [Sh, Ih, Rh] = simulate_absir(n_population,Iv0, T, infection_rate, recovery_rate,pod_size,gathering_chance)
% Simulate agent-based transmission model with pods and other cool stuff
%
% Inputs
%   n_population: Population size
%   Iv0 (column vector): Initial infection state
%   T (integer): Number of timesteps to simulate
%   infection_rate (float): Infection rate (probabilistic)
%   recovery_rate (float): Recovery rate (probabilistic)
%   pod_size (integer): Size of the pods
%   gathering_chance (float): Chance for agent to gather
% 
%
% Returns
%   Sh (matrix): Susceptible state history
%   Ih (matrix): Infected state history
%   Rh (matrix): Recovered state history

    % Setup
    Ih = zeros(n_population, T);  %Make infection history
    Ih(:, 1) = Iv0;                        %Record the initial state to infection history
    Rh = zeros(n_population, T); %Make recovery history
        
    %Define scrambler
    function [I] = scramble(I,order)
        for p=1:(length(I))
            out(p,1) = I(order(p));
        end
        I = out;
    end
    
    %Define unscrambler
    function [I] = unscramble(I,order)
        for p=1:(length(I))
            out(order(p),1) = I(p);
        end
        I = out;
    end

    % Construct an "action" helper function
    function [I, R,n_population,pod_size,gathering_chance] = action(I, R,n_population,pod_size,gathering_chance)

        %Choose a random set of people to remove
        remove = rand(1,n_population)>gathering_chance;
        
        %Store the data of those people
        store(:,1) = I(remove,:); store(:,2) = R(remove,:);

        %Remove the people from the data that goes into the simulation
        I(remove,:) = [];
        R(remove,:) = [];

        %scramble I & R
        order = randperm(length(I));
        I = scramble(I,order);
        R = scramble(R,order);

        %Make the social graph for the data that is going to be simulated
        M = pod_maker(size(I,1),pod_size);

        % Compute infection probabilities based on the social graph
        v_eff = infection_rate * M * I;
        v_pr = v_eff ./ (1 + v_eff);
        
        % Draw random values
        v_infect = rand(size(I,1), 1) <= v_pr;
        v_recover = rand(n_population, 1) <= recovery_rate;

        % Infect non-recovered individuals
        I = I | v_infect & (~R);
        new_infection = v_infect & (~R);

        %Unscramble I & R & new_infection
        I = unscramble(I,order);
        R = unscramble(R,order);
        new_infection = unscramble(new_infection, order);

        %merge people back together
        for g=1:n_population
            if remove(g) == 1
                newR(g,1) = store(1,2);
                newI(g,1) = store(1,1);
                store(1,:) = [];
                just_infected(g,1) = 0;
            else
                newR(g,1) = R(1);
                newI(g,1) = I(1);
                just_infected(g,1) = new_infection(1);
                R(1) = [];
                I(1) = [];
                new_infection(1) = [];
            end
        end
        I = newI;
        R = newR;

        % Recover infected individuals (added part that makes sure cant get infected and recover in same timestep)
        R = R | ((I & v_recover) & ~just_infected);
        I = I & (~R);

        %Assign outputs
        n_population = n_population;
        pod_size = pod_size;
        gathering_chance = gathering_chance;
    end
        
    % Run simulation
    for i = 2:T
        [Ih(:, i), Rh(:, i),n_population,pod_size,gathering_chance] = action(Ih(:, i-1), Rh(:, i-1),n_population,pod_size,gathering_chance);
    end
    
    % Compute susceptible history
    Sh = ones(n_population, T) - Ih - Rh;
end

