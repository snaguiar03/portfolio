## Project 1: How to find π?

The following project utilizes a Monte Carlo simulation, using varying levels of precision and graphically depicting random points generated within a unit square and differentiating between points that are inside and outside of the unit circle

Using the Monte Carlo method, the project estimates π via random point generation within a square with side length 2. The points inside in the circle are tracked and the ratio of the inside points to the total number of points allows us to approximate π. These calculations continue until the difference between the consecutive estimates of π fall below the precision level threshold. 

# Features: 

- Monte Carlo method: approximates π via random point generation and tracks the numner of points outside of circle
- Precision Control: specified levels expressed via significant figures (3)

# Usage: 
- Clone/download repo
- Use MATLAB to access folder with files
- Run precise_computer with desire sig figs as input

## Project 2: Mandelbrot Boundary Approximation

The project is divided into four parts: generating a fractal using iterative methods, approximating its boundary with a bisection algorithm, fitting a polynomial to the boundary, and calculating the boundary length with numerical integration. For the first part, the function fractal(c) calculates the number of iterations required for a given point c in the complex plane to diverge in accordance to z = z^2+c. When the magnitude of z exceeds 2 within 100 iterations, the points diverge. As a result, the function is applied to a grid of complex points spanning x ∈ [-2,1] and y ∈ [-2,1]. 

The bisection function approximates the boundary of the fractal by finding the point where an indicator function switches from negative (inside the fractal) to positive (outside the fractal). This method is applied along vertical slices of the complex plane, and the boundary points for each x-coordinate are found. The indicator function evaluates whether a point is inside or outside the fractal. In the polynomial fitting function, a 15th-order polynomial is fitted to the data using the boundary points from the bisection method. This process smooths the jagged boundary of the fractal for further analysis. The polynomial is fitted to a selected range of boundary points, discarding outliers on the left and right where the boundary flattens. Finally, the length of the polynomial boundary is computed using numerical integration. The approximate length of the fractal boundary comes out to 3.3061063. 

The polynomial approximation provides a smooth model of the boundary, though high-order polynomials can introduce artifacts outside the fitted range. The choice of polynomial order and range fitting was critical to avoid overfitting or misbehavior in boundary regions.
