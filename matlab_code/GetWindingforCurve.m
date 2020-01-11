function winding = GetWindingforCurve( curve ) 
%%%GETWINDINGFORCURVE get the winding number of a curve (presumablly in a plane)
%
% curve is a matrix of size [ dim , Nk ], where dim is the dimension of the ambient space, and Nk is
% the number of grid points for parametrizing the curve. This 1D curve lives in a dim-dimensional
% ambient space. 
%
% winding is an integer of the winding of the curve
%
% We define the winding number of the 1D curve in a dim-dimensional ambient space as follows.
% Suppose the curve is in a plane (we shall determine if this is the case later), then we shall
% first determine the (dim-2) vectors that are orthogonal to the curve. With that, we shall perform
% a coordinate transformation such that the curve lies within the x-y plane. Note that we didn't
% assume the plane of the curve cross the origin of the ambient space. After that, we can use the
% usual way to determine the winding of the curve. 

dim = size( curve , 1 ) ; 
Nk = size( curve , 2 ) ; 
dk = 2*pi/Nk ; 
krange = -pi : dk : (pi-dk) ; 

%% First take three points on the curve, and determines all the (dim-2) vectors that are orthorgonal to them
p1 = curve( : , 1 ) ; 
p2 = curve( : , floor( Nk/3 ) ) ; 
p3 = curve( : , 2 * floor( Nk/3 ) ) ; 
vec = [ p1 - p2 , p2 - p3 ] ; 
[Q,~] = qr( [ vec , zeros(dim,dim-2) ] , 0 ) ; 
vec_ortho = Q( : , 3 : end ) ; % The (dim-2) vectors that are orthogonal to the plane of the curve

% The other two basis vector, which span the plane of the curve are the following
vec(:,1) = vec(:,1) / norm( vec(:,1) ) ; 
vec(:,2) = vec(:,2) / norm( vec(:,2) ) ; 
% We pick the basis such that the first point is in the first quadrant.
if vec(:,1)' * p1 < 0
    vec(:,1) = -vec(:,1) ; 
end
if vec(:,2)' * p1 < 0
    vec(:,2) = -vec(:,2) ; 
end

newbasis = [ vec , vec_ortho ] ;  % The new basis

%% Now check that indeed the curve is in a plane
% % We start with point p1, and form vectors with all the other points, and check the orthogonality
% overlap = zeros( dim-2 , Nk-1 ) ; 
% for z = 1 : Nk-1 
%     overlap( : , z ) = vec_ortho' * ( p1 - curve( : , z+1 ) ) ; 
% end
% 
% if max( max( abs( overlap ) ) ) > 0.1 % 0.0001
%     winding = NaN ; 
%     figure ; 
%     plot( krange(2:end) , abs( overlap ) ) ; 
%     error( 'You are doomed, this curve is not on a plane' ) ; 
% end

%% If we pass the above test, then let's figure out the new representation of the curve in the new basis
h = zeros( 2 , Nk ) ; % The new representaton of the curve in the new basis
for z = 1 : Nk 
    h( : , z ) = newbasis( : , 1:2 )' * curve( : , z ) ; 
end
winding = -GetWinding( h(1,:) , h(2,:) ) ; 

end