function[X,Y] = EinthovenCoord(T1,T2)
% conver T to XY on Einthoven triangle

% Jan-2020  MA

%% internalParams
d120 = 2*pi*120/360;   % 120 degries
d30 = 2*pi*30/360;
s30 = sin(d30);
c120 = cos(d120);
s120 = sin(d120);

%% initialize
if ~all(size(T1)==size(T2))
    error('Matlab:Statistx:ImproperParam',...
        'T1 and T2 must have the same size')
end
Y = T1;
X = zeros(size(Y));
Sgn = ones(size(T1));
% Sgn(T1>0) = -1;
% Sgn(T2>0) = -1;
Sgn(T1.*T2<0) = -1;

%% compute  y = ax+b
a = s120/c120;  % slope
b = Sgn.*double(T2)/s30;    % intercept
% plot both both Positive or Both Negative
I = Sgn>0;
X(I) = (double(T1(I))+b(I))/a;
% convert differen signs
I = Sgn<=0;
X(I) = (double(T1(I))-b(I))/a;

return
