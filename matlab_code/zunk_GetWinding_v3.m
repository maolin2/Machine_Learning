function winding = GetWinding_v3( P )
%%%GETWINDING_V3 get the nontrivial winding of a non-trivial occupied manifold
%
% P is a matrix of [occ,occ,Nk] for ONE state in the occ manifold, where Nk is the number of k-point
% in the BZ

occ = size( P , 1 ) ;
Nk = size( P , 3 ) ;
% dk = 2*pi / Nk ;
% krange = -pi : dk : (pi-dk) ;

% By default, occ = 2^D, where the dimension of the Pauli matrices are of size [2^D,2^D], and the
% space of pauli matrices of dim 2D+1. First let's construct the Gamma matrices.
D = log(occ)/log(2) ;

Gamma = GetGammaMatrices( D ) ; % Gamma is a matrix of size [2^D , 2^D , 2*D+1]

d = zeros( 2*D+1 , Nk ) ; % The hamiltonian is H = d_i Gamma_i
H = zeros( 2^D , 2^D , Nk ) ;
for z = 1 : Nk
    for zz = 1 : 2*D+1
        d(zz,z) = trace( P(:,:,z) * Gamma( : , : , zz ) ) ;
        H( : , : , z ) = H( : , : , z ) + d( zz , z ) * Gamma( : , : , zz ) ;
    end
end

%%% Now we have assumed that the occupied manifold of P has nontrivial winidng. In other word, in
%%% the appropriate chosen basis, H should be represented as a 1D curve in the (2*D+1)-dim space of
%%% the Gamma matrices. Note that the current basis of Gamma need not be a good one. In other word,
%%% in the good basis, there is only two nonzero entries in the d-vector.


d2 = zeros( 2*D+1 , 2*D ) ; 
% There are 2*D vectors of length (2*D+1) need to be figured out. The 1st one would be one of the
% Gamma matrix associated with the 1D curve. The rest of the 2D-1 ones would be the one that is
% orthogonal to it. We shall figure out the second one Gamma matrix associated with the 1D curve, by
% multiplpying up these 2D matrix. We don't want to figure out this last matrix directly because in
% the QR process, we don't know if that will shuffle the order of the matrices. The only thing we
% are sure is that the first matrix will not be moved. 
for z = 1 : 2*D
    for zz = 1 : 2*D+1
        d2( zz , z ) = trace( H( : , : , z ) * Gamma( : , : , zz ) ) ;
    end
end

% Here we Gram-Schmidt orthonormalization the d2 vector. By default, d2(:,1) is not changed but
% normalized. Note, we know d2 has to be a real vector, thus we set it to be real here. If we not
% take only the real part, the qr decomposition may give complex entries, which will be wrong.
[d2, ~] = qr(real(d2),0);


Gamma2 = zeros( size(Gamma) ) ; % The good basis
temp = eye( 2^D ) ; % This would be the last Gamma2 matrix. 
for z = 1 : 2*D
    for zz = 1 : (2*D+1)
        Gamma2( : , : , z ) = Gamma2( : , : , z )  + d2( zz , z ) * Gamma( : , : , zz ) ;
    end
    temp = temp * Gamma2( : , : , z ) ; 
end
Gamma2( : , : , end ) = (-1i)^D * temp ; 

hx = zeros( 1 , Nk ) ; 
hy = zeros( 1 , Nk ) ; 
for z = 1 : Nk
    hx(z) = 1/(2^D) * trace( H(:,:,z) * Gamma2(:,:,1) ) ; 
    hy(z) = 1/(2^D) * trace( H(:,:,z) * Gamma2(:,:,end) ) ; 
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