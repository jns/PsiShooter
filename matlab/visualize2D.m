function visualize2D(varargin)
%function visualize2D(varargin)
%visualize2D([1:100],[1:100],peaks(100),0,[1:20],[1:20],peaks(20)+10,5);
%input data structure:
% varargin{1+3*n} = x basis vector ie: [x1,x2,x3]
% varargin{2+3*n} = y basis vector ie: {y1,y2,y3]
% varargin{3+3*n} = 2D data array ie: [d_11,d_12,d_13;d_21,d_22,d_23]
% level{4+3*n} = plot offset from zero energy (for the wavefunctions)
%where n is the number of items being plotted
%style will pass matlab plot options directly to the plot command.

if mod(nargin,4) ~= 0 && nargin ~= 0;
    fprintf(['Bad input arguments. There should be a multiple \n' ...
        'of four input arguments. There are ' num2str(nargin) '.\n']);
    return;
end

%plot the wavefunctions and potential together.
figure(1);
for n = 1:4:nargin
    if n == 1 %plot styles set up for the potential (more transparent)
        h(1)=surf(varargin{n},varargin{n+1},varargin{n+2});
        set(h(1),'facealpha',0.35);
        set(h(1),'edgealpha',0.35);
        hold on;
    else %plot styles set up for wavefunctions
        %h(n) is the handle for the graphic objects I am creating. They can
        %be accessed by get or set(h(n)) with whatever property you want to
        %mess with.
        h(n)=surf(varargin{n},varargin{n+1},varargin{n+2}+varargin{n+3});
        set(h(n),'facealpha',0.7);
        set(h(n),'edgealpha',0.7);
    end
end
hold off;
'end'% a good place to stick a breakpoint for playing with the plots.
