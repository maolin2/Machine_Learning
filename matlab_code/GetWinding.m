function winding = GetWinding( hx , hy )
%%% GETWINDING get the winding number of a curve parameterized by two vectors hx and hy
%
% hx and hy are vectors of size [1,Nk]
%
% Nk is the number of points in the 1D BZ

% dk = 2*pi/Nk ; 
% krange = -pi : dk : (pi-dk) ; 

U = hx + 1i * hy ; 

angleU = angle( U ) ; 

diff_angleU = [ angleU( 2 : end ) - angleU( 1 : end-1 ) , angleU(1) - angleU(end) ] ; 

epsilon = pi ; % This is the value beyond which we say there is a winding. This is picked randomly for now

winding = zeros( size( diff_angleU ) ) ; 
winding( abs(diff_angleU) > epsilon ) = sign( diff_angleU( abs(diff_angleU) > epsilon ) ) ; 

% figure ; 
% plot( 1 : length( diff_angleU ) , diff_angleU , 1 : length( diff_angleU ) , winding ) ; 

winding = sum( winding ) ; 

end