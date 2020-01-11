tic ; 

p0 = [ 1 0 ; 0 1 ];
p1 = [ 0 1 ; 1 0 ];
p2 = [ 0 -1i ; 1i 0 ];
p3 = [ 1 0 ; 0 -1 ];

G0 = kron( p3 , p0 ) ; 
G1 = -kron( p2 , p1 ) ;
G2 = -kron( p2 , p2 ) ;
G3 = -kron( p2 , p3 ) ;
G4 = kron( p1 , p0 ) ;

G5 = kron( p1 , p1 ) ;
G6 = kron( p1 , p2 ) ;

G52m = 1/2 * ( G5 - G2 ) ;
G52p = 1/2 * ( G5 + G2 ) ;
G16m = 1/2 * ( G1 - G6 ) ;
G16p = 1/2 * ( G1 + G6 ) ;

coef = rand(1,2) * 4 - 2 ; 
phirange = linspace( 0.01 , pi , 30 ) ; 
delrange = phirange ;
Nx = 20 ; 
Ny = 20 ; 

Hk = @(phi,kx,ky , del ) ( coef(1) + cos(kx) ) * G4 + sin(kx) * G3 ...
    + ( cos(phi) * ( coef(2) + cos(ky) ) + sin(phi) * sin(ky) ) * G52m ...
    + ( coef(2) + cos(ky) ) * G52p ...
    + ( sin(phi) * ( coef(2) + cos(ky) ) - cos(phi) * sin(ky) ) * G16p ...
    + sin(ky) * G16m + del * G0 ; 

Hx_amd = @(phi,kx , del ) ( coef(1) + cos(kx) ) * G4 + sin(kx) * G3 ...
    + ( cos(phi) * ( coef(2) ) ) * G52m ...
    + ( coef(2) ) * G52p ...
    + ( sin(phi) * ( coef(2) ) ) * G16p  + del * G0 ;


Hx_asub = @(phi) ( cos(phi) * ( 1/2 ) + sin(phi) * 1/(2i) ) * G52m ...
    + ( 1/2 ) * G52p ...
    + ( sin(phi) * ( 1/2 ) - cos(phi) * 1/(2i) ) * G16p ...
    + 1/(2i) * G16m ; 

Hx = @(kx,phi,del) blktridiag( Hx_amd(phi,kx,del) , Hx_asub(phi) , Hx_asub(phi)' , Nx ) ; 

Hy_amd = @(phi,ky,del) ( coef(1) ) * G4 ...
    + ( cos(phi) * ( coef(2) + cos(ky) ) + sin(phi) * sin(ky) ) * G52m ...
    + ( coef(2) + cos(ky) ) * G52p ...
    + ( sin(phi) * ( coef(2) + cos(ky) ) - cos(phi) * sin(ky) ) * G16p ...
    + sin(ky) * G16m + del * G0  ;

Hy_asub = @(phi) ( 1/2 ) * G4 + 1/(2i) * G3 ;
            
Hy = @(ky,phi,del) blktridiag( Hy_amd(phi,ky,del) , Hy_asub(phi) , Hy_asub(phi)' , Ny ) ; 

%%% x inside, y outside
Hobc_amd_amd = @(phi,del) ( coef(1) ) * G4 ...
    + ( cos(phi) * ( coef(2) ) ) * G52m ...
    + ( coef(2) ) * G52p ...
    + ( sin(phi) * ( coef(2) ) ) * G16p + del * G0 ; 

Hobc_amd = @(phi,del) blktridiag( Hobc_amd_amd(phi,del) , Hy_asub(phi) , Hy_asub(phi)' , Nx ) ; 

Hobc_asub = @(phi) blktridiag( Hx_asub(phi) , zeros(4) , zeros(4) , Nx ) ; 

Hobc = @(phi,del) blktridiag( Hobc_amd(phi,del) , Hobc_asub(phi) , Hobc_asub(phi)' , Ny ) ; 

% mx = 0 ; 
% my = 0 ; 
% gx = 1 ; 
% gy = 1 ; 
% Phi = pi ; 
% [Hk , Hx , Hy , Hobc , Hpbc] = GetQuadModels( mx , my , gx , gy , Phi , Nx , 1 ) ; 

EPx = zeros( 1, length(phirange) ) ; 
EPy = zeros( 1, length(phirange) ) ; 
Qcor = zeros( 1 , length(phirange) ) ; 

phi = pi/2 ; 
del = 0.001 ; 

parfor z = 1 : length(phirange) 
    
    disp( z ) ; 
    
    Hx0 = @(kx) Hx( kx , phirange(z) , del ) ; 
%     Hx0 = @(kx) Hx( kx , phi , delrange(z) ) ; 
    poldensx = pol2D( Hx0 , Nx , 4 , 100 , 'usual' ) ;
    EPx(z) = sum( poldensx( 1 : Nx/2 ) ) ; 
    
    Hy0 = @(ky) Hy( ky , phirange(z) , del ) ; 
%     Hy0 = @(ky) Hy( ky , phi , delrange(z) ) ; 
    poldensy = pol2D( Hy0 , Ny , 4 , 100 , 'usual' ) ;
    EPy(z) = sum( poldensy( 1 : Ny/2 ) ) ; 
    
    Hobc0 = Hobc( phirange(z) , del ) ; 
%     Hobc0 = Hobc( phi , delrange(z) ) ; 
    charge_dist = cor2D( Hobc0 , Nx , 4 ) ; 
    Qcor(z) = sum( sum( charge_dist( 1 : Nx/2 , 1 : Ny/2 ) ) ) ; 
    
end

EPx = mod( EPx , 1 ) ; 
EPy = mod( EPy , 1 ) ; 
Qcor = -mod( Qcor , 1 ) ; 

qxy = EPx + EPy - Qcor ;
qxy = mod( qxy , 1 ) ; 

figure ; 
plot( phirange , qxy , 'ro' , phirange , EPx , 'bx' , phirange , EPy , 'b*-' , phirange , Qcor , 'kx' ) ;
legend( 'qxy' , 'EPx' , 'EPy' , 'Qcor' ) ;
title( coef ) ; 
axis tight ; 
grid on

%%
% H0 = @(kx,ky) Hk(pi/2 , kx , ky , del ) ; 
% [ ~ , ~ , Wxy_occ , ~ , Wyx_occ, ~ ] = quad2D( H0 , 4 , 100 ) ;
% 
% figure ; 
% hold ; 
% plot( -abs(Wxy_occ) , 'ro' ) ; 
% plot( abs(Wyx_occ) , 'bx' ) ; 
% title( coef ) ; 

toc 