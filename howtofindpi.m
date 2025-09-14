%% Part 1: using monte carlo to approximate pi 

clear; clc; 

%randomize points for calculation
data_set = [10, 100, 1000, 5000, 10000, 100000, 500000, 1000000]; 

%exact value of pi
exact = pi;

%create arrays to store results

%monte carlo simulation

for idx = 1:length(data_set)
    N = data_set(idx); 

    tic; %begin timing simulation

    x = rand(N,1);
    y = rand(N,1);

    points_incircle = sum(x.^2+ y.^2 <= 1); %use period to square each element of x/y instead of matrix exponentiation
    
    %estimate the value of pi at this point
    estimate(idx) = 4 * points_incircle/N;


    %stop execution timer and record time
    execution_time(idx) = toc;

    %measure deviation from true value 
    deviation(idx) = abs(estimate(idx) - exact);

end

%Plots

% pi estimate
figure; 
subplot(2,1,1);
plot(points_incircle, estimate, 'b-o', 'LineWidth', 2);
hold on;
%using exact value as reference 
yline(exact, 'r--', 'LineWidth', 2);
xlabel('Total # of Points');
ylabel('Estimate')
title('Monte Carlo Simulation to Estimate Pi');
legend('Computed pi', 'True pi');
grid on; 

%plotting deviations
subplot(2,1,2);
plot(points_incircle, deviation, 'g-o', 'LineWidth',2);
xlabel('Number of Points');
ylabel('Deviation from true value');
title("Deviation of Estimate from True Value"); 
grid on; 

%computation cost vs precision Figure
figure; 
plot(points_incircle, deviation, 'g-o', 'LineWidth',2);
xlabel('Number of Points');
ylabel('Execution time (secs)');
title('Computational cost for Varying # of Points');
grid on; 


%% Part 2: using a while loop to specified precision level

% desired number of significant figures
sf = 3;

% tol based on significant figures
tol = 5 * 10^-(sf);  % Set tolerance to half the next significant figure

% Initialize
N = 10000; 
estimate = 0; 
delta_pi = Inf;  
prev_pi = 0;  
iterations = 0;  

% while loop to continue until the difference between successive estimates is within tol
while delta_pi > tol
    x = rand(N,1);
    y = rand(N,1);
    
    inside_circle = sum(x.^2 + y.^2 <= 1);
    
    estimate = 4 * inside_circle / N;
    
    delta_pi = abs(estimate - prev_pi);
    
    prev_pi = computed_pi;
    
    N = N + 1000;  
    iterations = iterations + 1;
end

fprintf('Computed value of π: %.10f\n', estimate);
fprintf('Number of iterations required: %d\n', iterations);
fprintf('Final number of points used: %d\n', N);

for idx = 1:length(N_values)
    % Repeat the pi computation for each N_value
    x = rand(N_values(idx), 1);
    y = rand(N_values(idx), 1);
    points_incircle = sum(x.^2 + y.^2 <= 1);
    pi_values(idx) = 4 * inside_circle / N_values(idx);
end

%evolution of computed π
figure;
plot(N_values, pi_values, 'b-o', 'LineWidth', 2);
yline(pi, 'r--', 'LineWidth', 2);  % Reference line for the true value of π
xlabel('Number of Points');
ylabel('Computed π');
title('Convergence of Monte Carlo Approximation of π');
grid on;
legend('Computed π', 'True π');


%% Part 3: Modifying part 2 

function computed_pi = precise_compute(sf)

    tolerance = 5 * 10^-(sf); 

    % Initialize variables
    N = 1000; 
    computed_pi = 0; 
    delta_pi = Inf;
    prev_pi = 0; 
    iteration_count = 0; 
    
    figure;
    hold on;
    axis equal;
    grid on;
    title('Monte Carlo Approximation of π');
    xlabel('x');
    ylabel('y');
    while delta_pi > tolerance
        x = rand(N,1) * 2 - 1;  % Random x in [-1, 1]
        y = rand(N,1) * 2 - 1;  % Random y in [-1, 1]
        
        inside_circle = (x.^2 + y.^2 <= 1);
        
        plot(x(inside_circle), y(inside_circle), 'b.');  % Points inside the circle in blue
        plot(x(~inside_circle), y(~inside_circle), 'r.');  % Points outside the circle in red
        drawnow;
        
        computed_pi = 4 * sum(inside_circle) / N;
        
        delta_pi = abs(computed_pi - prev_pi);
        
        prev_pi = computed_pi;
        
        N = N + 1000;  
        
        iteration_count = iteration_count + 1;
    end

    %final computed value of π w/ sig figs
    final_pi_str = sprintf(['Computed value of π to ', num2str(sig_figs), ' significant figures: %.', num2str(sig_figs), 'f\n'], computed_pi);
    fprintf(final_pi_str);
    
    % reference circle
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k-', 'LineWidth', 2);
    
    % Adding final pi computed value 
    text(-1.5, 1.2, ['\pi ≈ ', sprintf(['%.', num2str(sig_figs), 'f'], computed_pi)], 'FontSize', 14, 'Color', 'black');

    legend('Points inside circle', 'Points outside circle', 'Unit circle');
    
    return;
end
