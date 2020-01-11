% Here we study the winding of wannier bands for quadrupole models. 

% clearvars

tic;

n = 4;

Nkx = 100 ; 
Nky = 100 ; 
Nx = 10 ; % Not important for now
Ny = 10 ; % Not important for now
dkx = 2*pi/Nkx ; 
dky = 2*pi/Nky ; 
kxrange = -pi : dkx : (pi-dkx) ; 
kyrange = -pi : dky : (pi-dky) ; 


p0 = [ 1 0 ; 0 1 ];
p1 = [ 0 1 ; 1 0 ];
p2 = [ 0 -1i ; 1i 0 ];
p3 = [ 1 0 ; 0 -1 ];
G1 = -kron( p2 , p1 ) ; 
G2 = -kron( p2 , p2 ) ; 
G3 = -kron( p2 , p3 ) ; 
G4 = kron( p1 , p0 ) ; 
G5 = kron( p3 , p0 ) ; 
del = 0.00000001 ; 

%% Get Hamiltonian 
Nxc = 2 ; % The truncation of the fourier series along x
Nyc = 1 ; % The truncation of the fourier series along x
t1 = rand( 1 , Nxc+1 ) ; % The coefficients for the cos(kx) in G4
t2 = rand( 1 , Nxc+1 ) ; % The coefficients for the sin(kx) in G3
t3 = rand( 1 , Nyc+1 ) ; % The coefficients for the cos(ky) in G2
t4 = rand( 1 , Nyc+1 ) ; % The coefficients for the sin(ky) in G1

% t1 = [1/2,1] ; 
% t2 = [0,1] ; 
% t3 = [1/2,1] ; 
% t4 = [0,1] ; 

hx1 = @(kx) 0 ; 
hx2 = @(kx) 0 ; 
for z = 1 : Nxc+1 
    hx1 = @(kx) hx1(kx) + t1(z) * cos( (z-1) * kx ) ;
    hx2 = @(kx) hx2(kx) + t2(z) * sin( (z-1) * kx ) ;
end
hy1 = @(ky) 0 ; 
hy2 = @(ky) 0 ; 
for z = 1 : Nyc+1 
    hy1 = @(ky) hy1(ky) + t3(z) * cos( (z-1) * ky ) ;
    hy2 = @(ky) hy2(ky) + t4(z) * sin( (z-1) * ky ) ;
end

Hk = @(kx,ky) hx1(kx) * G4 + hx2(kx) * G3 + hy1(ky) * G2 + hy2(ky) * G1 + del * G5 ; 

%% Get Wannier states
% wx, wy are wannier states, which are of size [n,occ,Nkx,Nky]
[wx,wy] = GetWannierStates( Hk , Nkx , Nky ) ;

%% Calculate the NWL for checking
py_wx = zeros( 1 , Nkx ) ; % The nwl of wx(kx,ky) along ky, as a function of kx, whihc is a const
px_wy = zeros( 1 , Nky ) ; % The nwl of wy(kx,ky) along kx, as a function of ky, whihc is a const

parfor x = 1 : Nkx
    y0 = 1 ;
    loop = 1 ;
    for ii = 1 : Nky
        nexty0 = y0 + 1 ;
        if nexty0 == (Nky+1)
            nexty0 = 1 ;
        end
        loop = wx( : , 1 , x , nexty0 )' * wx( : , 1 , x , y0 ) * loop ;
        y0 = nexty0 ;
    end
    py_wx(x) = angle( eig( loop ) ) / (2*pi) ;
end

parfor y = 1 : Nky
    x0 = 1 ;
    loop = 1 ;
    for ii = 1 : Nkx
        nextx0 = x0 + 1 ;
        if nextx0 == (Nkx+1)
            nextx0 = 1 ;
        end
        loop = wy( : , 1 , nextx0 , y )' * wy( : , 1 , x0 , y ) * loop ;
        x0 = nextx0 ;
    end
    px_wy(y) = angle( eig( loop ) ) / (2*pi) ;
end

figure ; 
subplot(1,2,1) ; 
plot( kxrange , py_wx , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'NWL py wx(kx) along ky' ) ; 
xlabel('kx') ; 
subplot(1,2,2) ; 
plot( kyrange , px_wy , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'NWL px wy(ky) along kx' ) ; 
xlabel('ky') ; 

%% Get the smooth gauge, then calculate wx^\dagger i\partial_y wx
winding_wx = zeros( 1 , Nkx ) ; % The winding of wx(kx,ky) along ky, as a function of kx, whihc is a const
winding_wy = zeros( 1 , Nky ) ; % The winding of wy(kx,ky) along kx, as a function of ky, whihc is a const

Ay_wx = zeros( 1 , Nky ) ; 
for x = 1 : Nkx
    
    wx2 = GetParallelTransportGauge( permute( wx( : , 1 , x , : ) , [1,2,4,3] ) , 1 ) ;
    
    for y = 1 : Nky
        if y == 1
            Ay_wx(1) = 1i * wx2( : , 1 , 1 )' * ( wx2( : , 1 , 2 ) - wx2( : , 1 , Nky ) ) / (2*dky) ;
        elseif y == Nky
            Ay_wx(Nky) = 1i * wx2( : , 1 , Nky )' * ( wx2( : , 1 , 1 ) - wx2( : , 1 , Nky-1 ) ) / (2*dky) ;
        else
            Ay_wx(y) = 1i * wx2( : , 1 , y )' * ( wx2( : , 1 , y+1 ) - wx2( : , 1 , y-1 ) ) / (2*dky) ;
        end
    end
    
    winding_wx(x) = 2 * sum( Ay_wx ) / Nky ; 
end

Ax_wy = zeros( 1 , Nkx ) ; 
for y = 1 : Nky
    
    wy2 = GetParallelTransportGauge( wy( : , 1 , : , y ) , 1 ) ;
    
    for x = 1 : Nkx
        if x == 1
            Ax_wy(1) = 1i * wy2( : , 1 , 1 )' * ( wy2( : , 1 , 2 ) - wy2( : , 1 , Nkx ) ) / (2*dkx) ;
        elseif x == Nkx
            Ax_wy(Nkx) = 1i * wy2( : , 1 , Nkx )' * ( wy2( : , 1 , 1 ) - wy2( : , 1 , Nkx-1 ) ) / (2*dkx) ;
        else
            Ax_wy(x) = 1i * wy2( : , 1 , x )' * ( wy2( : , 1 , x+1 ) - wy2( : , 1 , x-1 ) ) / (2*dkx) ;
        end
    end
    
    winding_wy(y) = 2 * sum( Ax_wy ) / Nkx ; 
end


%%
% winding_wx = zeros( 1 , Nkx ) ; % The winding of wx(kx,ky) along ky, as a function of kx, whihc is a const
% winding_wy = zeros( 1 , Nky ) ; % The winding of wy(kx,ky) along kx, as a function of ky, whihc is a const
% 
% for z = 1 : Nkx
%     winding_wx(z) = GetWinding_v5( permute( wx( : , : , z , : ) , [1,2,4,3] ) ) ; 
% end
% for z = 1 : Nky
%     winding_wy(z) = GetWinding_v5( wy( : , : , : , z ) ) ; 
% end

%% Get the spin expectation evaluated for wx and wy respectively
% D = 2 ;
% Gamma = GetGammaMatrices( D ) ; % Gamma is a matrix of size [2^D , 2^D , 2*D+1]
% spinExp_wx = zeros( 2*D+1 , Nky , Nkx ) ; 
% spinExp_wy = zeros( 2*D+1 , Nkx , Nky ) ; 
% for z = 1 : Nkx
%     for zz = 1 : Nky
%         for zzz = 1 : 2*D+1
%             spinExp_wx( zzz , zz , z ) = trace( wx( : , : , z , zz )' * Gamma( : , : , zzz ) * wx( : , : , z , zz ) ) ; 
%         end
%     end
% end
% 
% for z = 1 : Nky
%     for zz = 1 : Nkx
%         for zzz = 1 : 2*D+1
%             spinExp_wy( zzz , zz , z ) = trace( wy( : , : , zz , z )' * Gamma( : , : , zzz ) * wy( : , : , zz , z ) ) ; 
%         end
%     end
% end
% winding_wx = zeros( 1 , Nkx ) ; % The winding of wx(kx,ky) along ky, as a function of kx, whihc is a const
% winding_wy = zeros( 1 , Nky ) ; % The winding of wy(kx,ky) along kx, as a function of ky, whihc is a const
% 
% for z = 1 : Nkx
%     winding_wx(z) = GetWindingforCurve( spinExp_wx( : , : , z ) ) ; 
% end
% for z = 1 : Nky
%     winding_wy(z) = GetWindingforCurve( spinExp_wy( : , : , z ) ) ; 
% end

%% Calculate the winding of corresponding 1DTI along x and y direction for comparison
windingx = GetWinding( hx1(kxrange) , hx2(kxrange) ) ;
windingy = GetWinding( hy1(kyrange) , hy2(kyrange) ) ;
disp( [ '[ windingy , windingx ] = ' , num2str([ windingy , windingx ]) ] ) ; 

%% Plot

figure ; 
subplot(1,2,1) ; 
plot( kxrange , real(winding_wx) , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'Winding of wx(kx) along ky' ) ; 
xlabel('kx') ; 


subplot(1,2,2) ; 
plot( kyrange , real(winding_wy) , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'Winding of wy(ky) along kx' ) ; 
xlabel('ky') ; 
disp( t1 ) ; 
disp( t2 ) ; 
disp( t3 ) ; 
disp( t4 ) ; 


%% Get the smooth gauge, then calculate wx^\dagger i\partial_x wx, The unconventional NWL
winding_wx_y = zeros( 1 , Nky ) ; 
winding_wy_x = zeros( 1 , Nkx ) ; 

Ax_wx = zeros( 1 , Nkx ) ; 
for y = 1 : Nky
    
    wx2 = GetParallelTransportGauge( wx( : , 1 , : , y ) , 1 ) ;
    
    for x = 1 : Nkx
        if x == 1
            Ax_wx(1) = 1i * wx2( : , 1 , 1 )' * ( wx2( : , 1 , 2 ) - wx2( : , 1 , Nkx ) ) / (2*dkx) ;
        elseif x == Nkx
            Ax_wx(Nkx) = 1i * wx2( : , 1 , Nkx )' * ( wx2( : , 1 , 1 ) - wx2( : , 1 , Nkx-1 ) ) / (2*dkx) ;
        else
            Ax_wx(x) = 1i * wx2( : , 1 , x )' * ( wx2( : , 1 , x+1 ) - wx2( : , 1 , x-1 ) ) / (2*dkx) ;
        end
    end
    
    winding_wx_y(y) = 2 * sum( Ax_wx ) / Nkx ; 
end

Ay_wy = zeros( 1 , Nky ) ; 
for x = 1 : Nkx
    
    wy2 = GetParallelTransportGauge( permute( wy( : , 1 , x , : ) , [1,2,4,3] ) , 1 ) ;
    
    for y = 1 : Nky
        if y == 1
            Ay_wy(1) = 1i * wy2( : , 1 , 1 )' * ( wy2( : , 1 , 2 ) - wy2( : , 1 , Nky ) ) / (2*dky) ;
        elseif y == Nky
            Ay_wy(Nky) = 1i * wy2( : , 1 , Nky )' * ( wy2( : , 1 , 1 ) - wy2( : , 1 , Nky-1 ) ) / (2*dky) ;
        else
            Ay_wy(y) = 1i * wy2( : , 1 , y )' * ( wy2( : , 1 , y+1 ) - wy2( : , 1 , y-1 ) ) / (2*dky) ;
        end
    end
    
    winding_wy_x(x) = 2 * sum( Ay_wy ) / Nky ; 
end

figure ; 
subplot(1,2,1) ; 
plot( kyrange , real(winding_wx_y) , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'Winding of wx(ky) along kx' ) ; 
xlabel('kx') ; 


subplot(1,2,2) ; 
plot( kxrange , real(winding_wy_x) , 'o' ) ; 
axis tight ; 
grid on ; 
title( 'Winding of wy(kx) along ky' ) ; 
xlabel('ky') ; 
disp( t1 ) ; 
disp( t2 ) ; 
disp( t3 ) ; 
disp( t4 ) ;

toc ; 