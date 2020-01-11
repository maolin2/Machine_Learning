function [BPx,APx,EPx,BPy,APy,EPy,Qcor] = GetBoundaryQuantities_v2(mx , my , gx , gy , Phi , Nx , Ny , del,Nk,model)
%%%GETBOUNDARYQUANTITIES get the following boundaries quantities for a given quadrupole models
%
% The differnce compared to v1 is that we allow Nx not the same as Ny
%
% BPx, bulk pol along x, a vector of [Ndel,1]
% APx, actual edge pol along x, a vector of [Ndel,2]
% EPx, edge pol along x, a vector of [Ndel,1]
% BPy, bulk pol along y, a vector of [Ndel,1]
% APy, actual edge pol along y, a vector of [Ndel,2]
% EPy, edge pol along y, a vector of [Ndel,1]
% Qcor, corner mode, a vector of [Ndel,1]


[~ , Hx , Hy , Hobc , ~] = GetQuadModels_v2( mx , my , gx , gy , Phi , Nx , Ny , model ) ; 
% Hk(kx,ky,del) is a functional handel of 4x4 matrix 
% Hx(kx,del) is a functional handle of [4*N,4*N] matrix 
% Hy(ky,del) is a functional handle of [4*N,4*N] matrix 
% Hobc(del) is a functional handle of [4*N^2,4*N^2] matrix

Ndel = length(del) ; 
n = 4 ; 

BPxdens = zeros( Ndel , Ny ) ;
BPydens = zeros( Ndel , Nx ) ;
APxdens = zeros( Ndel , 2,Ny ) ;
APydens = zeros( Ndel , 2,Nx ) ;
EPxdens = zeros( Ndel , Ny ) ;
EPydens = zeros( Ndel , Nx ) ;
Qcor = zeros( Ndel , 1 ) ;
for z = 1 : Ndel
    disp(z) ;
    H0 = @(kx) Hx( kx , del(z) ) ;
    [ BPxdens(z,:) , APxdens(z,:,:) , EPxdens(z,:) ] = pol2D_v2( H0 , Ny , n , Nk , 'usual' ) ;
    
    H0 = @(ky) Hy( ky , del(z) ) ;
    [ BPydens(z,:) , APydens(z,:,:) , EPydens(z,:) ] = pol2D_v2( H0 , Nx , n , Nk , 'usual' ) ;
    
    H0 = Hobc( del(z) ) ;
    charge_dist = cor2D_v2( H0 , Nx , Ny , n ) ;
    Qcor(z) = -sum(sum( charge_dist( 1:Nx/2 , 1:Ny/2 )-2) ) ;
end

BPx = sum( BPxdens( : , 1:end/2 ) , 2 ) ; 
BPy = sum( BPydens( : , 1:end/2 ) , 2 ) ; 
APx = sum( APxdens( : , : , 1:end/2 ) , 3 ) ; 
APy = sum( APydens( : , : , 1:end/2 ) , 3 ) ; 
EPx = sum( EPxdens( : , 1:end/2 ) , 2 ) ; 
EPy = sum( EPydens( : , 1:end/2 ) , 2 ) ;  

for z = 1 : Ndel
    if EPx(z) > 0
        EPx(z) = EPx(z) - 1 ;
    end
    if EPy(z) > 0
        EPy(z) = EPy(z) - 1 ;
    end
    if Qcor(z) > 0
        Qcor(z) = Qcor(z) - 1 ;
    end
end
% quad_def = -1-( Edge_pol_x + Edge_pol_y - cordist.' ) ;


end