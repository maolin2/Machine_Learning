% Here we try to understand the higher winding of 1DTI in the basis of px and py
% H = \sum_n ( t1(n) cos(n*k) ) * p1 + \sum_n ( t2(n) sin(n*k) ) * p2 
%
% We have three ways to check the winding number.
% 1. Calculate hx, hy on the complex plane. Check the winding as we vary k.
%
% 2. Calculate using wilson loop. But this only works for winding number = 1. In particular, if
% winding = 1, we have 1i * log(W)/(2*pi) = 1/2, and winding = 0 we have that be 0.
%
% 3. Calculate using w = 1i/(2*pi) \int dk U(k)^* ( \partial_k U(k) ) 
%
% We note that using the endpoint to calculate it is wrong, due to the singularity

% clearvars 

p0 = [ 1 0 ; 0 1 ];
p1 = [ 0 1 ; 1 0 ];
p2 = [ 0 -1i ; 1i 0 ];
p3 = [ 1 0 ; 0 -1 ];

N = 1 ; % The highest order in the fourier series

Nk = 1000 ; 
dk = 2*pi/Nk ; 
krange = -pi : dk : (pi-dk) ; 

t1 = rand( 1 , N+1 ) ; % The coefficients for the cos in px
t2 = rand( 1 , N+1 ) ; % The coefficients for the sin in px
t3 = rand( 1 , N+1 ) ; % The coefficients for the cos in py
t4 = rand( 1 , N+1 ) ; % The coefficients for the sin in py

%% Construct the model 
hx = @(k) 0 ; 
hy = @(k) 0 ; 
for n = 0 : N 
    hx = @(k) hx(k) + t1(n+1) * cos( n * k ) + t2(n+1) * sin( n * k ) ;
    hy = @(k) hy(k) + t3(n+1) * cos( n * k ) + t4(n+1) * sin( n * k ) ;
end

% hx = @(k) cos( 1 * k ) ; 
% hy = @(k) sin( 1 * k ) ; 

%% Forget about other methods, simply smooth the gauge, then take the differivative
H = @(k) hx(k) * p1 + hy(k) * p2 ; 

eigvec = zeros( 2 , 1 , Nk ) ;
A = zeros( 1 , Nk ) ; 
for z = 1 : Nk
    [V,D] = eig( H( krange(z) ) ) ; 
    [D,I] = sort( diag(D) ) ; 
    V = V(:,I) ; 
    eigvec( : , : , z ) = V( : , 1 ) ; 
end
eigvec2 = GetParallelTransportGauge( eigvec  , 1 ) ;
% eigvec2 = eigvec ; 

for z = 1 : Nk
    if z == 1
        A(1) = 1i * eigvec2( : , : , 1 )' * ( eigvec2( : , : , 2 ) - eigvec2( : , : , Nk ) ) / (2*dk) ;
    elseif z == Nk
        A(Nk) = 1i * eigvec2( : , : , Nk )' * ( eigvec2( : , : , 1 ) - eigvec2( : , : , Nk-1 ) ) / (2*dk) ;
    else
        A(z) = 1i * eigvec2( : , : , z )' * ( eigvec2( : , : , z+1 ) - eigvec2( : , : , z-1 ) ) / (2*dk) ;
    end
end

winding = 2 * sum( A ) / Nk ; 

%% Get the Winding, using the naive basis H = hx px + hy py
% winding = GetWinding( hx(krange) , hy(krange) ) ;

%% Get the winding from the spin expectation evaluated with eigenstates
% H = zeros( 2 , 2 , Nk ) ;
% for z = 1 : Nk
%     H( : , : , z ) = hx(krange(z)) * p1 + hy(krange(z)) * p2 ;
% end
% V = zeros( 2 , 1 , Nk ) ;
% for z = 1 : Nk
%     [Vtemp , D] = eig( H( : , : , z ) ) ;
%     [~,I] = sort( D ) ;
%     Vtemp = Vtemp( : , I ) ; 
%     V( : , : , z ) = Vtemp( : , 1 ) ; 
% end
% 
% D = 1 ;
% Gamma = GetGammaMatrices( D ) ; % Gamma is a matrix of size [2^D , 2^D , 2*D+1]
% spinExp = zeros( 2*D+1 , Nk ) ; 
% for z = 1 : Nk
%     for zz = 1 : 2*D+1
%         spinExp( zz , z ) = trace( V( : , : , z )' * Gamma( : , : , zz ) * V( : , : , z ) ) ; 
%     end
% end
% 
% winding3 = GetWindingforCurve(spinExp) ; 

%% Plot the winding 
figure ; 
scatter( hx(krange) , hy(krange) , 5 ) ; 
grid on ; 
axis( [-2,2,-2,2] ) ; 
title( num2str( [ real(winding) ] ) ) ; 
