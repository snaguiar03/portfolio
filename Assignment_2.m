% Mandelbrot Boundary Approximation
clear; close all; clc;

%% parameters
max_iter = 100;        
N_plot   = 200;        
n_boundary = 1000;     
tol_bisect = 1e-6;

%% Part 1: Mandelbrot visualization
x_vis = linspace(-2, 1, N_plot);
y_vis = linspace(-1, 1, N_plot);
[X, Y] = meshgrid(x_vis, y_vis);
C = X + 1i * Y;

IT = zeros(N_plot, N_plot);
for ii = 1:N_plot
    for jj = 1:N_plot
        IT(ii, jj) = fractal(C(ii,jj), max_iter);
    end
end

figure;
imagesc(x_vis, y_vis, IT); axis xy; colorbar;
title('Mandelbrot iteration counts'); xlabel('x'); ylabel('y');

%% Part 2: Extract boundary points
xvals = linspace(-2, 1, n_boundary);
boundary_y = nan(size(xvals));

for k = 1:length(xvals)
    xx = xvals(k);
    fn = indicator_fn_at_x(xx, max_iter);
    Ys = linspace(-1, 1, 401);
    vals = arrayfun(fn, Ys);
    idx = find(diff(sign(vals)) ~= 0, 1);

    if ~isempty(idx)
        s = Ys(idx); e = Ys(idx+1);
        boundary_y(k) = bisection(fn, s, e, tol_bisect);
    end
end

valid = ~isnan(boundary_y);
xvals = xvals(valid);
boundary_y = boundary_y(valid);

figure;
plot(xvals, boundary_y, '.');
title('Extracted boundary points'); xlabel('x'); ylabel('y');

%% Part 3: Clean up boundary for polynomial fitting
% remove big jumps in y
dy = [0 diff(boundary_y)];
mask = abs(dy) < 0.2; % threshold
xvals = xvals(mask);
boundary_y = boundary_y(mask);

% keep largest continuous chunk
gapIdx = find(diff(xvals) > 0.01);
edges = [0, gapIdx, numel(xvals)];

segLengths = diff(edges);
[~, largestIdx] = max(segLengths);

segStart = edges(largestIdx)+1;
segEnd   = edges(largestIdx+1);

%largest contiguous segment
xfit = xvals(segStart:segEnd);
yfit = boundary_y(segStart:segEnd);

figure;
plot(xfit, yfit, 'bo');
title('Cleaned boundary points for fitting'); xlabel('x'); ylabel('y');

%% polynomial fit
order = 15;
p = polyfit(xfit, yfit, order);
xpoly = linspace(min(xfit), max(xfit), 1000);
ypoly = polyval(p, xpoly);
figure;
plot(xfit, yfit, 'bo', 'MarkerSize', 4); hold on;
plot(xpoly, ypoly, 'r-', 'LineWidth', 1.2);
legend('Boundary points','15th-order fit','Location','Best');
title('Fractal boundary with polynomial fit');
xlabel('x'); ylabel('y');

%% Part 4: Length of polynomial boundary
s = min(xfit);
e = max(xfit);
len = poly_len(p, s, e);

fprintf('Approximate length of the fractal boundary (poly fit): %.6f\n', len);

%% Local functions
function it = fractal(c, max_iter)
    z = 0;
    it = max_iter;
    for n = 1:max_iter
        z = z^2 + c;
        if abs(z) > 2
            it = n;
            return;
        end
    end
end

function m = bisection(fn_f, s, e, tol)
    if nargin < 4, tol = 1e-6; end
    fs = fn_f(s); fe = fn_f(e);
    if fs * fe > 0
        m = NaN; return;
    end
    while (e - s) > tol
        m = (s + e) / 2;
        fm = fn_f(m);
        if fm == 0, break; end
        if sign(fm) == sign(fe)
            e = m; fe = fm;
        else
            s = m; fs = fm;
        end
    end
    m = (s + e) / 2;
end

function fn = indicator_fn_at_x(x, max_iter)
    fn = @(y) double(fractal(x + 1i*y, max_iter) < max_iter) * 2 - 1;
end

function l = poly_len(p, s, e)
    dp = polyder(p);
    ds = @(x) sqrt(1 + (polyval(dp, x)).^2);
    l = integral(ds, s, e, 'RelTol', 1e-8, 'AbsTol', 1e-10);
end
