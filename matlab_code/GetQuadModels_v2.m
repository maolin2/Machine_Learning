function [Hk , Hx , Hy , Hobc , Hpbc] = GetQuadModels_v2( mx , my , gx , gy , Phi , Nx , Ny , model )
%%% GETQUADDMODELS returns four Hamiltonians of quadrupole models, in different geometries
%
% % The differnce compared to v1 is that we allow Nx not the same as Ny
%
% Nx and Ny are the system size
% Hk(kx,ky,del) is a functional handel of 4x4 matrix
% Hx(kx,del) is a functional handle of [4*N,4*N] matrix
% Hy(ky,del) is a functional handle of [4*N,4*N] matrix
% Hobc(del) is a functional handle of [4*N^2,4*N^2] matrix
% Hpbc(del) is a functional handle of [4*N^2,4*N^2] matrix
%
%
% C4T, Phi not applicable
%         Hk = @(kx,ky,del)  sin(kx) * G3 + sin(ky) * G1 + del * G0 ...
%             + (-2 + mx + cos(kx) + cos(ky) ) * G4  + my * ( cos(kx) - cos(ky) ) * G2 ;

%%
n = 4;
if model ~= 5
    
    Hk = @(kx,ky,del) [ del , 0 , my + gy * exp(1i*ky) , mx + gx * exp(1i*kx) ; ...
        0 , del , mx + gx * exp(-1i*kx) , exp(1i*Phi) * ( my + gy * exp(-1i*ky) ) ; ...
        my + gy * exp(-1i*ky) , mx + gx * exp(1i*kx) , -del , 0 ; ...
        mx + gx * exp(-1i*kx) , exp(-1i*Phi) * ( my + gy * exp(1i*ky) ) , 0 , -del ] ;
    
    Hyamd = @(ky,del) [ del , 0 , my + gy * exp(1i*ky) , mx ; ...
        0 , del , mx  , exp(1i*Phi) * ( my + gy * exp(-1i*ky) ) ; ...
        my + gy * exp(-1i*ky) , mx , -del , 0 ; ...
        mx , exp(-1i*Phi) * ( my + gy * exp(1i*ky) ) , 0 , -del ] ;
    Hyasub = [ 0 , 0 , 0 , 0 ; ...
        0 , 0 , gx , 0 ; ...
        0 , 0 , 0 , 0 ; ...
        gx , 0 , 0 , 0 ] ;
    Hy = @(ky,del) blktridiag( Hyamd(ky,del) , Hyasub , Hyasub' , Nx ) ;
    
    Hxamd = @(kx,del) [ del , 0 , my , mx + gx * exp(1i*kx) ; ...
        0 , del , mx + gx * exp(-1i*kx) , exp(1i*Phi) * ( my ) ; ...
        my , mx + gx * exp(1i*kx) , -del , 0 ; ...
        mx + gx * exp(-1i*kx) , exp(-1i*Phi) * ( my ) , 0 , -del ] ;
    Hxasub = [ 0 , 0 , 0 , 0 ; ...
        0 , 0 , 0 , exp(1i*Phi) * gy ; ...
        gy , 0 , 0 , 0 ; ...
        0 , 0 , 0 , 0 ] ;
    Hx = @(kx,del) blktridiag( Hxamd(kx,del) , Hxasub , Hxasub' , Ny ) ;
    
    %%% x inside, y outside
    Hobc_amd_amd = @(del) [ del , 0 , my , mx ; ...
        0 , del , mx , exp(1i*Phi) * (my) ; ...
        my , mx , -del , 0 ; ...
        mx , exp(-1i*Phi) * (my) , 0 , -del ] ;
    
    Hobc_amd_asub = [ 0 , 0 , 0 , 0 ; ...
        0 , 0 , gx , 0 ; ...
        0 , 0 , 0 , 0 ; ...
        gx , 0 , 0 , 0 ] ;
    Hobc_amd = @(del) blktridiag( Hobc_amd_amd(del) , Hobc_amd_asub , Hobc_amd_asub' , Nx ) ;
    
    Hobc_asub_amd = [ 0 , 0 , 0 , 0 ; ...
        0 , 0 , 0 , exp(1i*Phi) * gy ; ...
        gy , 0 , 0 , 0 ; ...
        0 , 0 , 0 , 0 ] ;
    Hobc_asub = blktridiag( Hobc_asub_amd , zeros(n) , zeros(n) , Nx ) ;
    
    Hobc = @(del) blktridiag( Hobc_amd(del) , Hobc_asub , Hobc_asub' , Ny ) ;
    
    % impose PBC along x
    amd_pbc = zeros( 4*Nx , 4*Nx ) ;
    amd_pbc( 1 : 4 , end-3:end ) = Hobc_amd_asub ;
    amd_pbc( end-3:end , 1 : 4 ) = (Hobc_amd_asub)' ; 
    Hpbc_amd = @(del) Hobc_amd(del) + amd_pbc ;
    Hpbc = @(del) blktridiag( Hpbc_amd(del) , Hobc_asub , Hobc_asub' , Ny ) ; 
    
    % impose PBC along y 
    pbc = zeros( 4 * Nx * Ny , 4 * Nx * Ny ) ; 
    pbc( 1 : (4*Nx) , ( end - 4*Nx + 1 ) : end ) =  Hobc_asub ;
    pbc( ( end - 4*Nx + 1 ) : end , 1 : (4*Nx) ) =  Hobc_asub' ;
    Hpbc = @(del) Hpbc(del) + pbc ; 
    
elseif model==5
    [p0,p1,p2,p3] = GetPauliMatrices();
    G0 = kron( p3 , p0 );
    G1 = kron( -p2 , p1 );
    G2 = kron( -p2 , p2 );
    G3 = kron( -p2 , p3 );
    G4 = kron( p1 , p0 );
    
    Hk = @(kx,ky,del)  sin(kx) * G3 + sin(ky) * G1 + del * G0 ...
        + (-2 + mx + cos(kx) + cos(ky) ) * G4  + my * ( cos(kx) - cos(ky) ) * G2 ;
    
    Hyamd = @(ky,del) sin(ky) * G1 + del * G0 + (-2 + mx + cos(ky) ) * G4  + my * ( - cos(ky) ) * G2 ;
    Hyasub = G3/(2i) + G4/2 + (my) * G2/2 ; 
    Hy = @(ky,del) blktridiag( Hyamd(ky,del) , Hyasub , Hyasub' , Nx ) ;
    
    Hxamd = @(kx,del) sin(kx) * G3 + del * G0 + (-2 + mx + cos(kx) ) * G4  + my * ( cos(kx) ) * G2 ;
    Hxasub = G1/(2i) + G4/2 + (-my) * G2/2 ; 
    Hx = @(kx,del) blktridiag( Hxamd(kx,del) , Hxasub , Hxasub' , Ny ) ;
    
     %%% x inside, y outside
    Hobc_amd_amd = @(del) del * G0 + (-2 + mx ) * G4 ;
    
    Hobc_amd_asub = G3/(2i) + G4/2 + my * G2/2 ; 
    Hobc_amd = @(del) blktridiag( Hobc_amd_amd(del) , Hobc_amd_asub , Hobc_amd_asub' , Nx ) ;
    
    Hobc_asub_amd = G1/(2i) + G4/2 + (-my) * G2/2 ;
    Hobc_asub = blktridiag( Hobc_asub_amd , zeros(n) , zeros(n) , Nx ) ;
    
    Hobc = @(del) blktridiag( Hobc_amd(del) , Hobc_asub , Hobc_asub' , Ny ) ;
%     crap = 0.01 * rand(size(Hobc(rand))) ; 
%     crap = 1/2 * ( crap + crap' ) ; 
%     Hobc = @(del) Hobc(del) + crap ; 
    
    % impose PBC along x
    amd_pbc = zeros( 4*Nx , 4*Nx ) ;
    amd_pbc( 1 : 4 , end-3:end ) = Hobc_amd_asub ;
    amd_pbc( end-3:end , 1 : 4 ) = (Hobc_amd_asub)' ; 
    Hpbc_amd = @(del) Hobc_amd(del) + amd_pbc ;
    Hpbc = @(del) blktridiag( Hpbc_amd(del) , Hobc_asub , Hobc_asub' , Ny ) ; 
    
    % impose PBC along y 
    pbc = zeros( 4 * Nx * Ny , 4 * Nx * Ny ) ; 
    pbc( 1 : (4*Nx) , ( end - 4*Nx + 1 ) : end ) =  Hobc_asub ;
    pbc( ( end - 4*Nx + 1 ) : end , 1 : (4*Nx) ) =  Hobc_asub' ;
    Hpbc = @(del) Hpbc(del) + pbc ; 
end


end
