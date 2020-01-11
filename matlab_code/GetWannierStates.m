function [wx,wy] = GetWannierStates( Hk , Nkx , Nky ) 
%%%GETWANNIERSTATES get the Wannier states for a given Hk(kx,ky) for BZ of size Nkx , Nky
%
% Hk(kx,ky) is a functional handle
%
% Nkx, Nky is the dim of the 2D BZ
%
% wx, wy are wannier states, which are of size [n,occ,Nkx,Nky] where n is the size of the
% Hamiltonand and occ is the number of filled bands

dkx = 2*pi/Nkx ;
dky = 2*pi/Nky ;
n = size( Hk(rand,rand) , 1 ) ; 
occ = n / 2 ; 

kxrange = -pi : dkx : (pi-dkx) ; 
kyrange = -pi : dky : (pi-dky) ; 

eigvec = zeros( 4 , 2 , Nkx , Nky ) ; 

for x = 1 : Nkx
    for y = 1 : Nky
        [V,D] = eig( full( Hk( kxrange(x) , kyrange(y) ) ) ) ; 
        [~,I] = sort( diag(D) ) ; 
        V = V(:,I) ; 
        eigvec( : , : , x , y ) = V( : , 1 : 2 ) ; 
    end
end

%% Calculate ( the exponential of ) Wannier Hamiltonian

eigvecWx = zeros( occ , occ , Nkx , Nky ) ;
parfor y = 1 : Nky
    for x = 1 : Nkx
        x0 = x ;
        loop = 1 ;
        for ii = 1 : Nkx
            nextx0 = x0 + 1 ;
            if nextx0 == (Nkx+1)
                nextx0 = 1 ;
            end
            loop = eigvec( : , : , nextx0 , y )' * eigvec( : , : , x0 , y ) * loop ;
            
%             temp = eigvec( : , : , nextx0 , y )' * eigvec( : , : , x0 , y ) ; 
%             [U,~,V] = svd( temp ) ; 
%             temp = U*V' ; 
%             loop = temp * loop ;

            x0 = nextx0 ;
        end
        [V,D1] = eig( loop );
        D1 = (angle(diag(D1)))/(2*pi) ; % Make sure D is within [0,1]
        [~,I] = sort( real( D1 ) ) ;
        eigvecWx(:,:,x,y) = V(:,I) ;
    end
end

eigvecWy = zeros( occ , occ , Nkx , Nky ) ;
parfor x = 1 : Nkx
    for y = 1 : Nky
        y0 = y ;
        loop = 1 ;
        for ii = 1 : Nky
            nexty0 = y0 + 1 ;
            if nexty0 == (Nky+1)
                nexty0 = 1 ;
            end
            loop = eigvec( : , : , x , nexty0 )' * eigvec( : , : , x , y0 ) * loop ;
            y0 = nexty0 ;
        end
        [V,D2] = eig( loop );
        D2 = (angle(diag(D2)))/(2*pi) ; % Make sure D is within [0,1]
        [~,I] = sort( real( D2 ) ) ;
        eigvecWy(:,:,x,y) = V(:,I) ;
    end
end

%% Calculate Wannier bands
wx = zeros( n , occ , Nkx , Nky ); % Linear combination of Bloch bands.
wy = zeros( n , occ , Nkx , Nky ); % Linear combination of Bloch bands.
for x = 1:Nkx
    for y = 1:Nky
        wx(:,:,x,y) = eigvec(:,:,x,y) * eigvecWx(:,:,x,y);
        wy(:,:,x,y) = eigvec(:,:,x,y) * eigvecWy(:,:,x,y);
    end
end

end