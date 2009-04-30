function visualize1D(varargin)
%function visualize1D(varargin)
%ex: visualize1D([1 2 3 4 5],[1 2 3 2 1],'red');
%input data structure:
% varargin{1+3*n} = X-axis array
% barargin{2+3*n} = Y-axis array
% style{3+3*n} = display style
%where n is the number of items being plotted
%style will pass matlab plot options directly to the plot command.

if mod(nargin,3) ~= 0 && nargin ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of three input arguments. There are ' num2str(nargin) '.\n']);
    return;
end

for n = 1:3:nargin
    X(:,ceil(n/3)) = varargin{n};
    Y(:,ceil(n/3)) = varargin{n+1};
    color{ceil(n/3)} = varargin{n+2};
end

figure(1);
for n = 1:1:nargin/3
    plot(X(:,n),Y(:,n),color{n})
    hold on;
end
hold off