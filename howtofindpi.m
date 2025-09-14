%% Part 1: using monte carlo to approximate pi 
clear; clc; close all;

data_set = [10, 100, 1000, 5000, 10000, 100000, 500000, 1000000]; 

exact = pi;
estimate = zeros(size(data_set));
execution_time = zeros(size(data_set));
deviation = zeros(size(data_set));

% monte carlo simulation
for idx = 1:length(data_set)
    N = data_set(idx); 

    tic; % begin timing simulation

    x = rand(N,1);
    y = rand(N,1);

    points_incircle = sum(x.^2 + y.^2 <= 1);
    
    % estimate the value of pi
    estimate(idx) = 4 * points_incircle / N;

    % stop execution timer and record time
    execution_time(idx) = toc;

    % measure deviation from true value 
    deviation(idx) = abs(estimate(idx) - exact);
end

% (a) π estimate vs N, with deviation subplot
figure; 
subplot(2,1,1);
plot(data_set, estimate, 'b-o', 'LineWidth', 2);
hold on;
yline(exact, 'r--', 'LineWidth', 2);
xlabel('Total # of Points');
ylabel('Estimate');
title('Monte Carlo Simulation to Estimate Pi');
legend('Computed π', 'True π');
grid on;

subplot(2,1,2);
plot(data_set, deviation, 'g-o', 'LineWidth', 2);
xlabel('Number of Points');
ylabel('Deviation from true value');
title('Deviation of Estimate from True Value');
grid on;

% (b) computational cost vs N
figure; 
plot(data_set, execution_time, 'm-o', 'LineWidth', 2);
xlabel('Number of Points');
ylabel('Execution time (secs)');
title('Computational Cost for Varying # of Points');
grid on;


%% Part 2: using a while loop to specified precision level

% desired number of significant figures
sf = 3;

tol = 5 * 10^(-sf);

N = 10000;  
estimate = 0; 
delta_pi = Inf;  
prev_pi = 0;  
iterations = 0;  

% arrays to track convergence
N_values = [];
pi_values = [];

% while loop until tolerance is met
while delta_pi > tol
    x = rand(N,1);
    y = rand(N,1);
    
    inside_circle = sum(x.^2 + y.^2 <= 1);
    
    estimate = 4 * inside_circle / N;
    
    delta_pi = abs(estimate - prev_pi);
    
    prev_pi = estimate;
    
    % record values for plotting
    N_values(end+1) = N; %#ok<*SAGROW>
    pi_values(end+1) = estimate;
    
    % increment N
    N = N + 1000;  
    iterations = iterations + 1;
end

% print results
fprintf('Computed value of π: %.10f\n', estimate);
fprintf('Number of iterations required: %d\n', iterations);
fprintf('Final number of points used: %d\n', N);

% plot convergence
figure;
plot(N_values, pi_values, 'b-o', 'LineWidth', 2);
yline(pi, 'r--', 'LineWidth', 2);  % reference line for true π
xlabel('Number of Points');
ylabel('Computed π');
title('Convergence of Monte Carlo Approximation of π');
grid on;
legend('Computed π', 'True π');


%% Part 3: Modifying part 2 

function computed_pi = precise_compute(sf)    
    tolerance = 5 * 10^(-sf);
    N = 1000; 
    computed_pi = 0; 
    delta_pi = Inf;
    prev_pi = 0; 
    iteration_count = 0; 
    
    figure;
    hold on;
    axis equal;
    grid on;
    title('Monte Carlo Approximation of \pi');
    xlabel('x');
    ylabel('y');
    
    while delta_pi > tolerance
        % generate random points in [-1,1]x[-1,1]
        x = rand(N,1) * 2 - 1;  
        y = rand(N,1) * 2 - 1;  
        
        % check which are inside circle
        inside_circle = (x.^2 + y.^2 <= 1);
        
        % plot points
        plot(x(inside_circle), y(inside_circle), 'b.');
        plot(x(~inside_circle), y(~inside_circle), 'r.');
        drawnow;
        
        % compute estimate
        computed_pi = 4 * sum(inside_circle) / N;
        
        % update stopping condition
        delta_pi = abs(computed_pi - prev_pi);
        prev_pi = computed_pi;
        
        % increment sample size
        N = N + 1000;  
        iteration_count = iteration_count + 1;
    end
    
    % print result to command window
    final_pi_str = sprintf(['Computed value of \x03C0 to %d significant figures: %.', ...
                             num2str(sf), 'f\n'], sf, computed_pi);
    fprintf(final_pi_str);
    
    % reference circle
    theta = linspace(0, 2*pi, 500);
    plot(cos(theta), sin(theta), 'k-', 'LineWidth', 2);
    
    % annotate figure with final pi
    text(-1.5, 1.2, ['\pi ≈ ', sprintf(['%.', num2str(sf), 'f'], computed_pi)], ...
        'FontSize', 14, 'Color', 'black');
    
    legend('Points inside circle', 'Points outside circle', 'Unit circle');
    
end

precise_compute(3);
