
function integ = integral(x,y)

%function integral, x, y
%
%;+
%; ROUTINE:        INTEGRAL
%;
%; USEAGE:         RESULT = INTEGRAL( X, Y )
%;
%; PURPOSE:        Integrate tabulated data using Simpson's rule
%;                 with 3-point Lagragian interpolation. Data may be
%;                 regularly sampled in X, or irregularly sampled.
%;
%; INPUT:
%;   X             Vector of x axis points.
%;                 (Elements must be unique and monotonically increasing)
%;
%;   Y             Vector of corresponding Y axis points.
%;
%; KEYWORD_INPUT:  None.
%;
%; OUTPUT:         Result of integration.
%;
%; EXAMPLE:
%;;       Example 1:
%;;       Define 11 x-values on the closed interval [0.0 , 0.8].
%;        X = [ 0.0, .12, .22, .32, .36, .40, .44, .54, .64, .70, .80 ]
%;
%;;       Define 11 f-values corresponding to x(i).
%;        F = [ 0.200000, 1.30973, 1.30524, 1.74339, 2.07490, 2.45600, $
%;               2.84299,  3.50730, 3.18194, 2.36302, 0.231964 ]
%;
%;;       Compute the integral.
%;        RESULT = INTEGRAL( X, F )
%;
%;;       In this example, the f-values are generated from a known function,
%;;       (f = .2 + 25*x - 200*x^2 + 675*x^3 - 900*x^4 + 400*x^5)
%;
%;;       The Multiple Application Trapezoid Method yields;  result = 1.5648
%;;       The Multiple Application Simpson's Method yields;  result = 1.6036
%;;              IDL User Library INT_TABULATED.PRO yields;  result = 1.6232
%;;                                    INTEGRAL.PRO yields;  result = 1.6274
%;;         The Exact Solution (4 decimal accuracy) yields;  result = 1.6405
%;
%; AUTHOR:         Liam Gumley, CIMSS/SSEC (liam.gumley@ssec.wisc.edu)
%;                 Based on a FORTRAN-77 version by Paul van Delst, CIMSS/SSEC
%;                 22-DEC-95
%;
%; REVISIONS:      None.
%;-
%
%; check that x and y have same number of elements
%
%n = n_elements( x )
n=length(x);
%if n ne n_elements( y ) then message, 'X and Y vectors must be same size'
%     
%; check that all points in X are unique, so we don't divide by zero
%
%x0 = double( x( 0 : n - 2 ) )
%x1 = double( x( 1 : n - 1 ) )
%loc = where( 1.0 + x1 - x0 le 1.0, count )
%if count ge 1 then begin
%  message, 'Non-unique X values were found at element numbers', /continue
%  print, strcompress( loc + 1, /remove_all )
%  message, ' '
%endif
%
%; make interpolation point arrays
%
%x0 = double( x( 0 : n - 3 ) )
%x1 = double( x( 1 : n - 2 ) )
%x2 = double( x( 2 : n - 1 ) )
%y0 = double( y( 0 : n - 3 ) )
%y1 = double( y( 1 : n - 2 ) )
%y2 = double( y( 2 : n - 1 ) )

x0 = x(1:n-2);
x1 = x(2:n-1);
x2 = x(3:n-0);
y0 = y(1:n-2);
y1 = y(2:n-1);
y2 = y(3:n-0);

%
%; compute interpolation delta and midpoint arrays
%
dx = x1 - x0;
xmid = 0.5.*( x1 + x0 );
%
%; compute 3 point lagrange interpolation
%
l0 = ((xmid-x1)./(x0-x1)).*((xmid-x2)./(x0-x2));
l1 = ((xmid-x0)./(x1-x0)).*((xmid-x2)./(x1-x2));
l2 = ((xmid-x0)./(x2-x0)).*((xmid-x1)./(x2-x1));
ymid = y0.*l0 + y1.*l1 + y2.*l2;

%
%; compute integral sum
%
integ = sum(1.0./6.0.*dx.*(y0+4.0.*ymid+y1));
%
%; handle last 3 points similarly
%
x0 = x( n - 2 );
x1 = x( n - 1 );
x2 = x( n - 0 );
y0 = y( n - 2 );
y1 = y( n - 1 );
y2 = y( n - 0 );
dx = x2 - x1;
xmid = 0.5.*( x2 + x1 );
l0 = ( ( xmid - x1 )./( x0 -x1 ) ).*( ( xmid - x2 )./( x0 - x2 ) );
l1 = ( ( xmid - x0 )./( x1 -x0 ) ).*( ( xmid - x2 )./( x1 - x2 ) );
l2 = ( ( xmid - x0 )./( x2 -x0 ) ).*( ( xmid - x1 )./( x2 - x1 ) );
ymid = y0.*l0 + y1.*l1 + y2.*l2;
integ = integ +  1.0/6.0.*dx.*( y1 + 4.0.*ymid + y2 );
%
%return, sum
%
%end




















