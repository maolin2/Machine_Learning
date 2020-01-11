clearvars

tic;
% clearvars ;

n = 4;
occ = n/2 ;
Nwocc = occ/2 ;
Nk = 30 ; % 30
dk = 2*pi/Nk ;
krange = -pi : dk : (pi-dk) ;

Ndel = 20 ;
del = linspace( -1 , -0.01 , Ndel ) ;
N = 20 ;

model = 3 ; % 1 mirror , 2 C4 flux , 3 Chiral , 4 Fake quad
%% Define models

if model == 1
    %%% mirror
    % Conclusion: Stripe mirror is Done.
    gx = 1.4 ;
    gy = 1. ;
    mx = 0.501 ;
    my = 0.001 ;
    Phi = pi ; % Fixed
elseif model == 2
    %%% C4 flux
    gx = 1 ;
    gy = gx ; % Fixed
    mx = 0. ;
    my = mx ; % Fixed
    Phi = pi/2 ;
    
elseif model == 3
    %%% Chiral model
    gx = 1.4 ;
    gy = 1. ;
    mx = 1.101 ;
    my = 0.501 ;
    Phi = pi/2 ;
elseif model == 4
    %%% Fake Quad
    gx = 1.5 ;
    gy = 1. ;
    mx = 1.1 ;
    my = 0. ;
    Phi = 0 ;  % Fixed
end

[Hk , ~ , ~ , ~,~] = GetQuadModels( mx , my , gx , gy , Phi , N , model ) ; 

[BPx,APx,EPx,BPy,APy,EPy,Qcor] = GetBoundaryQuantities( mx , my , gx , gy , Phi , N , del , Nk , model ) ; 

quad_def = ( EPx + EPy - Qcor ) ;

%%% Check for the decoupled cases
% Conclusion: It is correct. 
% a = sqrt( gx^2 + gy^2 + del.^2 - 2 * gx * gy * cos(Phi/2) ) ;
% b = sqrt( gx^2 + gy^2 + del.^2 + 2 * gx * gy * cos(Phi/2) ) ;
% z1 = ( gx - gy * exp(1i*Phi/2) ) ./ sqrt( 2 .* ( a + del ) .* a ) ;
% z3 = ( gx + gy * exp(1i*Phi/2) ) ./ sqrt( 2 .* ( b + del ) .* b ) ;
% z5 = sqrt( ( 1 -abs(z1).^2 ) .* ( 1 -abs(z3).^2 ) ) - conj(z1) .* z3 ;
% z6 = sqrt( ( 1 -abs(z1).^2 ) .* ( 1 -abs(z3).^2 ) ) + conj(z1) .* z3 ;
% U = z5.*conj(z6) + z6.*conj(z5) ; 
% figure ; 
% plot( del , quad_def , 'ro' , ...
%     del , ( -(1-abs(z5))./abs(z5) .* U/4 ) + ( -(1-abs(z6))./abs(z6) .* U/4 ) - ( 1/2 + del/4 .* (1./a+1./b) ) ) ; 
% axis tight ; 
% grid on ; 
% legend( 'quad def' , 'decoupled values' ) ; 

dk = 2*pi/Nk ;
krange2 = -pi : dk : (pi-dk) ;
[ qxyk , qyxk , Wyx_k , Wxy_k , Wx , Wy ] = quad2D_decoupled( Hk , n , Nk , del ) ; 
quad_yx_sum = permute( sum( sum( qyxk , 1 ) , 2 ) / (Nk^2) , [ 3 , 1 , 2 ] ) ; 
quad_xy_sum = permute( sum( sum( qxyk , 1 ) , 2 ) / (Nk^2) , [ 3 , 1 , 2 ] ) ; 

%% Get nested wilson loop 
% Nk2 = Nk ; 
% dk = 2*pi/Nk2 ;
% krange2 = -pi : dk : (pi-dk) ;
% Wx = zeros( Ndel , occ , Nk2 ) ; % as a function of ky
% Wy = zeros( Ndel , occ , Nk2 ) ; % as a function of kx
% Wxy_k = zeros( Ndel , occ , Nk2 ) ; % y first then x
% Wyx_k = zeros( Ndel , occ , Nk2 ) ; % x first then y
% Wxy = zeros( Ndel , occ ) ; % y first then x
% Wyx = zeros( Ndel , occ ) ; % x first then y
% 
% for z = 1 : Ndel
%     disp(z) ;
%     H = @(kx,ky) Hk( kx,ky , del(z) ) ;
%     [ Wx(z,:,:) , Wy(z,:,:) , Wxy_k( z , 1:Nwocc , : ) , Wxy_k( z , Nwocc+1:occ , : ) , ...
%         Wyx_k( z , 1:Nwocc , : ) , Wyx_k( z , Nwocc+1:occ , : ) ] = quad2D( H , n , Nk2 ) ; % Wxy is y first, then x
%     
%     for i = 1 : occ
%         for ii = 1 : Nk2
%             if Wxy_k(z,i,ii)>0
%                 Wxy_k(z,i,ii) = Wxy_k(z,i,ii)-1 ;
%             end
%             if Wyx_k(z,i,ii)>0
%                 Wyx_k(z,i,ii) = Wyx_k(z,i,ii)-1 ;
%             end
%         end
%     end
%     Wyx( z , : ) = sum( Wyx_k( z , : , : ) , 3 ) / Nk2 ;
%     Wxy( z , : ) = sum( Wxy_k( z , : , : ) , 3 ) / Nk2 ;
% end
% 
% %% Get eigvalWyWx + eigvalWy * eigvalWyWx AND eigvalWxWy + eigvalWx * eigvalWxWy
% quad_yx_sum = zeros( 1 , Ndel ) ;
% quad_xy_sum = zeros( 1 , Ndel ) ;
% for z = 1 : Ndel
%     for x = 1 : Nk2
%         for y = 1 : Nk2
%             quad_yx_sum( z ) = quad_yx_sum( z ) + Wyx_k( z , 2 , x ) ...
%                 + Wy( z , 1 , x ) * Wxy_k( z , 1 , y ) ...
%                 + Wy( z , 2 , x ) * Wxy_k( z , 2 , y ) ;
%             quad_xy_sum( z ) = quad_xy_sum( z ) + Wxy_k( z , 2 , y ) ...
%                 + Wx( z , 1 , x ) * Wyx_k( z , 1 , y ) ...
%                 + Wx( z , 2 , x ) * Wyx_k( z , 2 , y ) ;
%         end
%     end
%     quad_yx_sum(z) = quad_yx_sum(z) / Nk2^2 ;
%     quad_xy_sum(z) = quad_xy_sum(z) / Nk2^2 ;
% end
% 
% %%% Average eigvalWx(ky) and eigvalWy(kx) first 
% quad_yx_avg = zeros( 1 , Ndel ) ;
% quad_xy_avg = zeros( 1 , Ndel ) ;
% for z = 1 : Ndel
%     temp = 0 ; 
%     for x = 1 : Nk2
%         temp = temp + Wy(z,:,x) ; 
%     end
%     temp = temp / Nk2 ; 
%     quad_yx_avg(z) = Wyx(z,2) + temp(1) * Wxy(z,1) + temp(2) * Wxy(z,2) ; 
%     
%     temp = 0 ; 
%     for y = 1 : Nk2
%         temp = temp + Wx(z,:,y) ; 
%     end
%     temp = temp / Nk2 ; 
%     quad_xy_avg(z) = Wxy(z,2) + temp(1) * Wyx(z,1) + temp(2) * Wyx(z,2) ; 
% end

%% Plot

[X,Y] = meshgrid(krange2,del) ;

%%% Plot the dispersion of Wannier orb and nested wilson loop
figure ;
subplot(2,3,1) ;
hold ;
mesh( X , Y , permute( Wyx_k(:,2,:) , [1,3,2] ) ) ;
mesh( X , Y , permute( Wyx_k(:,1,:) , [1,3,2] ) ) ;
title( 'eigvalWyWx occ and eigvalWyWx unocc' ) ; 
xlabel('kx' ) ; 
ylabel('pump para') ;
view(90,0)
grid on
axis( [ krange2(1) , krange2(end) , del(1) , del(end) , -1 , 0 ] ) ;

subplot(2,3,2) ;
hold ;
mesh( X , Y , permute( Wxy_k(:,2,:) , [1,3,2] ) ) ;
mesh( X , Y , permute( Wxy_k(:,1,:) , [1,3,2] ) ) ;
title( 'eigvalWxWy occ and eigvalWxWy unocc' ) ; 
xlabel('ky' ) ; 
ylabel('pump para') ;
view(90,0)
grid on
axis( [ krange2(1) , krange2(end) , del(1) , del(end) , -1 , 0 ] ) ;

subplot(2,3,3) ;
hold ;
mesh( X , Y , permute( Wx(:,2,:) , [1,3,2] ) ) ;
mesh( X , Y , permute( Wx(:,1,:) , [1,3,2] ) ) ;
title( 'eigvalWx occ and eigvalWx unocc' ) ; 
xlabel('ky' ) ; 
ylabel('pump para') ;
view(90,0)
grid on
axis( [ krange2(1) , krange2(end) , del(1) , del(end) , -1/2 , 1/2 ] ) ;

subplot(2,3,4) ;
hold ;
mesh( X , Y , permute( Wy(:,2,:) , [1,3,2] ) ) ;
mesh( X , Y , permute( Wy(:,1,:) , [1,3,2] ) ) ;
title( 'eigvalWy occ and eigvalWy unocc' ) ; 
xlabel('kx' ) ; 
ylabel('pump para') ;
view(90,0)
grid on
axis( [ krange2(1) , krange2(end) , del(1) , del(end) , -1/2 , 1/2 ] ) ;

%%% Plot the distriubtion of Wannier Hamiltonian and Bulk dipole in cylinder geometry
% subplot(2,3,5) ;
% hold ;
% mesh( X , Y , permute( polx_bulk(1,:,:) , [ 3 , 2 , 1 ] ) ) ;
% % mesh( X , Y , permute( polx_bulk(2,:,:) , [ 3 , 2 , 1 ] ) ) ;
% xlabel('kx') ;
% ylabel('pump para') ;
% title( 'dispersion of pol x bulk' ) ;
% axis tight ;
% grid on
% view(90,0) ;
% 
% subplot(2,3,6) ;
% hold ;
% mesh( X , Y , permute( poly_bulk(1,:,:) , [ 3 , 2 , 1 ] ) ) ;
% % mesh( X , Y , permute( poly_bulk(2,:,:) , [ 3 , 2 , 1 ] ) ) ;
% xlabel('kx') ;
% ylabel('pump para') ;
% title( 'dispersion of pol y bulk' ) ;
% axis tight ;
% grid on
% view(90,0) ;

%% Plot the quad moment
figure ;
% plot( del , -1-quad_def , 'k' , del , quad_yx_sum , 'r*' , del , quad_xy_sum , 'bo' , del , quad_yx_avg , 'rv' , del , quad_xy_avg , 'b^' ) ;
% legend( 'quad def' , 'quad yx sum' , 'quad xy sum' , 'quad yx avg' , 'quad xy avg' ) ;
plot( del , -1-quad_def , 'k' , del , quad_yx_sum , 'r*' , del , quad_xy_sum , 'bo' ) ;
legend( 'quad def' , 'quad yx sum' , 'quad xy sum' ) ;
title( [ 'Everything, [mx , my,gx,gy,Phi,N,Nk] = ' , num2str( [mx , my,gx,gy,Phi,N,Nk] ) ] ) ;
grid on ;

toc

%% Get quad by summing Actual edge polarization with bulk polarization
%%% The two leg ladder Ham
% Hx_bulk = @(kx,del) [ -del , ( mx + gx*exp(1i*kx) ) , gy , 0 ; ...
%     ( mx + gx*exp(-1i*kx) ) , del , 0 , exp(1i*Phi) * gy ; ...
%     gy , 0 , del , mx + gx * exp(1i*kx) ; ...
%     0 , exp(-1i*Phi) * gy , mx + gx * exp(-1i*kx) , -del ] ; 
% 
% Hy_bulk = @(ky,del) [ -del , exp(-1i*Phi) * ( my + gy*exp(1i*ky) ) , gx , 0 ; ...
%     exp(1i*Phi) * ( my + gy*exp(-1i*ky) ) , del , 0 , gx ; ...
%     gx , 0 , del , my + gy * exp(1i*ky) ; ...
%     0 , gx , my + gy * exp(-1i*ky) , -del ] ;
% 
% eigvec = zeros( 4 , 2 , Nk ) ;
% WannierHx = zeros( 2 , 2 , Nk , Ndel ) ; % Wannier Hamiltonian along x
% WannierHy = zeros( 2 , 2 , Nk , Ndel ) ; % Wannier Hamiltonian along y
% eigvalWannierHx = zeros( 2 , Nk , Ndel ) ; % eigval of Wannier Ham along x
% eigvalWannierHy = zeros( 2 , Nk , Ndel ) ; % eigval of Wannier Ham along y
% polx_bulk = zeros( 2 , Nk , Ndel ) ; % The polarization ariss from the 2-leg ladder
% poly_bulk = zeros( 2 , Nk , Ndel ) ; % The polarization ariss from the 2-leg ladder
% for z = 1 : Ndel
%     for x = 1 : Nk
%         H20 = Hx_bulk( krange(x) , del(z) ) ;
%         [V,D] = eig( H20 ) ;
%         [D,I] = sort( diag(D) ) ;
%         V = V(:,I) ;
%         eigvec( : , : , x ) = V( : , 1:2 ) ;
%         %         eigvec( : , : , x ) = V( : , 1:2 ) * U(krange(x)) ;
%     end
%     for x = 1 : Nk
%         x0 = x ;
%         loop = 1 ;
%         for ii = 1 : Nk
%             nextx0 = x0 + 1 ;
%             if nextx0 == (Nk+1)
%                 nextx0 = 1 ;
%             end
%             loop = eigvec( : , : , nextx0 )' * eigvec( : , : , x0 ) * loop ;
%             x0 = nextx0 ;
%         end
%         WannierHx( : , : , x , z ) = (-1i) * logm( loop ) / (2*pi) ;
%         %         disp( WHx( : , : , x , z ) - WHx( : , : , x , z )' ) ;
%         [V,D] = eig( WannierHx( : , : , x , z ) );
%         [D,I] = sort( diag(D) ) ;
%         eigvalWannierHx( : , x , z ) = angle(D) ;
%         temp1 = eigvec( 1:2 , : , x ) * WannierHx( : , : , x , z ) * eigvec( 1:2 , : , x )' ;
%         temp2 = eigvec( 3:4 , : , x ) * WannierHx( : , : , x , z ) * eigvec( 3:4 , : , x )' ;
%         polx_bulk(2,x,z) = trace( temp1 ) ; % choice that consistent with eigvalWxWy
%         polx_bulk(1,x,z) = trace( temp2 ) ;
%     end
%     
%     for y = 1 : Nk
%         H20 = Hy_bulk( krange(y) , del(z) ) ;
%         [V,D] = eig( H20 ) ;
%         [D,I] = sort( diag(D) ) ;
%         V = V(:,I) ;
%         eigvec( : , : , y ) = V( : , 1:2 ) ;
%         %         eigvec( : , : , x ) = V( : , 1:2 ) * U(krange(x)) ;
%     end
%     for y = 1 : Nk
%         y0 = y ;
%         loop = 1 ;
%         for ii = 1 : Nk
%             nexty0 = y0 + 1 ;
%             if nexty0 == (Nk+1)
%                 nexty0 = 1 ;
%             end
%             loop = eigvec( : , : , nexty0 )' * eigvec( : , : , y0 ) * loop ;
%             y0 = nexty0 ;
%         end
%         WannierHy( : , : , y , z ) = (-1i) * logm( loop ) / (2*pi) ;
%         %         disp( WHy( : , : , y , z ) - WHy( : , : , y , z )' ) ;
%         [V,D] = eig( WannierHy( : , : , y , z ) );
%         [D,I] = sort( diag(D) ) ;
%         eigvalWannierHy( : , y , z ) = angle(D) ;
%         temp1 = eigvec( 1:2 , : , y ) * WannierHy( : , : , y , z ) * eigvec( 1:2 , : , y )' ;
%         temp2 = eigvec( 3:4 , : , y ) * WannierHy( : , : , y , z ) * eigvec( 3:4 , : , y )' ;
%         poly_bulk(2,y,z) = trace( temp1 ) ; % choice that consistent with eigvalWyWx
%         poly_bulk(1,y,z) = trace( temp2 ) ;
%     end
%     
% end
% polx_bulk = real(polx_bulk) ;
% poly_bulk = real(poly_bulk) ;
% 
% %%% We simply average the bluk dipole here.
% quad3yx = Wyx + permute( sum(poly_bulk,2) , [ 3 , 1 , 2 ] )  / Nk ;
% quad3xy = Wxy + permute( sum(polx_bulk,2) , [ 3 , 1 , 2 ] )  / Nk ;
% 
% % quad3yx = zeros( 2 , Ndel ) ; 
% % quad3xy = zeros( 2 , Ndel ) ; 
% % quad3yx(1,:) = eigvalWyWx(2,:) + sum( permute( poly_bulk(1,:,:) , [ 1 , 3 , 2 ] ) , 3 ) / Nk ;
% % quad3yx(2,:) = eigvalWyWx(1,:) + sum( permute( poly_bulk(2,:,:) , [ 1 , 3 , 2 ] ) , 3 ) / Nk ;
% % quad3xy(1,:) = eigvalWxWy(2,:) + sum( permute( polx_bulk(1,:,:) , [ 1 , 3 , 2 ] ) , 3 ) / Nk ;
% % quad3xy(2,:) = eigvalWxWy(1,:) + sum( permute( polx_bulk(2,:,:) , [ 1 , 3 , 2 ] ) , 3 ) / Nk ;
