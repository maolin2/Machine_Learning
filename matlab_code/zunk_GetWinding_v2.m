function winding = zunk_GetWinding_v2( H ) 
%%% GETWINDING_v2 get the winding number of an effective 1D Hamiltonian H(k) 
% This is a genearlized version of GetWinding, whic only works in a speciifc basis where H(k) = hx
% px + hy py
%
% H is a matrix of size [2,2,Nk] for the 1D BZ.


% p0 = [ 1 0 ; 0 1 ];
p1 = [ 0 1 ; 1 0 ];
p2 = [ 0 -1i ; 1i 0 ];
p3 = [ 1 0 ; 0 -1 ];

Nk = size(H,3) ; 

% dk = 2*pi/Nk ; 
% krange = -pi : dk : (pi-dk) ; 


%% First figure ONE basis
% Note that the basis is not unique, think of planar circle in the 3D sigma matrices space. The
% basis to write H(k) = hx sx + hy sy is not unique, where s_{x,y} are the effecive sigma matrices.
% But we can figure out one in the following way. We realize first that the circle is traced by
% evolving k from -pi to pi. Thus we can identifying sx as the vector corresponidng to H(pi) in the
% sigma space, which is a const matrix. Next, we randomly pick another matrix H(0) which has to be
% NOT parallel to H(pi) if we imagine them as vectors. But we note that H(0) need not be orthogonal
% to H(pi) either. Thus it cannot be another pauli matrix. However, with the two vectors H(pi) and
% H(0) we can figureout the vector that is perpendicular to the planar circle, which corresponds to
% sz. To determine sz, it its simply the cross product of the vectors of H(0) and H(pi). After we
% have sz, sy = [sz,sx]/(2i).

% a0 = trace( H(pi) * p0 ) ; 
a1 = trace( H(:,:,end) * p1 ) ; 
a2 = trace( H(:,:,end) * p2 ) ; 
a3 = trace( H(:,:,end) * p3 ) ; 
a = [a1,a2,a3] ; 

% The first effective Pauli matrix, which is represnted as the unit vector
% (a1,a2,a3)/sqrt(a1^2+a2^2+a3^2) in the pauli matrices space
s1 = ( a(1) * p1 + a(2) * p2 + a(3) * p3 ) /  norm(a) ; 


% Another vector in the pauli matrix space, which is NOT parallel to (a1,a2,a3)
% b0 = trace( H(0) * p0 ) ; 
b1 = trace( H(:,:,end/2) * p1 ) ; 
b2 = trace( H(:,:,end/2) * p2 ) ; 
b3 = trace( H(:,:,end/2) * p3 ) ; 

c = cross( [a1,a2,a3] , [b1,b2,b3] ) ; 

% The first effective Pauli matrix, which is represnted as the unit vector
% (c1,c2,c3)/sqrt(c1^2+c2^2+c3^2) in the pauli matrices space
s3 = ( c(1) * p1 + c(2) * p2 + c(3) * p3 ) / norm(c) ; 

% The second effective Pauli matrix
s2 = ( s3 * s1 - s1 * s3 ) / 2i ; 

hx = zeros( 1 , Nk ) ; 
hy = zeros( 1 , Nk ) ; 
for z = 1 : Nk
    hx(z) = 1/2 * trace( H(:,:,z) * s1 ) ; 
    hy(z) = 1/2 * trace( H(:,:,z) * s2 ) ; 
end

U = hx + 1i * hy ; 

angleU = angle( U ) ; 

diff_angleU = angleU( 2 : end ) - angleU( 1 : end-1 ) ; 

epsilon = pi ; % This is the value beyond which we say there is a winding. This is picked randomly for now

winding = zeros( size( diff_angleU ) ) ; 
winding( abs(diff_angleU) > epsilon ) = sign( diff_angleU( abs(diff_angleU) > epsilon ) ) ; 

% figure ; 
% plot( 1 : length( diff_angleU ) , diff_angleU , 1 : length( diff_angleU ) , winding ) ; 

winding = sum( winding ) ; 

end