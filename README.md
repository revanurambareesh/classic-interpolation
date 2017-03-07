# Interpolation

1. Cubic Spline Interpolation
2. Neville's Polynomial Interpolation
3. Linear curvefit
4. Polynomial 2nd order
5. Polynomial 3rd order
6. Polynomial 4th order

### Demo

Results:

![Console image](https://raw.githubusercontent.com/revanurambareesh/classic-interpolation/master/Resources/console1.png)

Depiction of f(x):

![Console image BGI](https://raw.githubusercontent.com/revanurambareesh/classic-interpolation/master/Resources/console2.png)


## ERROR LIST
- -9997 : Point not in range.
- -9996 : Number of points less than 4 (or greater than MAXSET)
- -9995 : Solution not converging in Newton Raphson's method
- -9994 : Set of equations in not diagonally dominant (or 0 element in matrix)
- -9993 : Linear backfit undefined for constant function. f(x) = c

##### Note

Folder 'PlotGraphs' contains different files created using Code::Blocks, among which Plotter.cpb must be chosen. The folder 'Resources' contains the WinBGI for displaying the plots generated.
The file PLOTTER.cpb must be opened in Code::Blocks.

The function plotCurve(int n){...} is where the actual Graphics is implemented. To generate a specific curve, lines could be commented out where putpixel(...) functions are called.
