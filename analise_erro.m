%% tcl link colors
%% $ns color 1 red
%% $ns color 2 green
%% $ns color 3 blue
%% $ns color 4 yellow
%% $ns color 5 magenta
%% $ns color 6 brown

red = [ 255 0 0 ];
green = [ 0 255 0 ];
blue = [ 0 0 255 ];
yellow = [ 255 255 0 ];
magenta = [ 255 0 255 ];
brown = [ 165 42 42 ];

color0 = red/255;
color1 = green/255;
color2 = blue/255;
color6 = yellow/255;
color7 = magenta/255;
color8 = brown/255;


FigHandle = figure;
set(FigHandle, 'Position', [0, 0, 640, 480]);
max_matrix_size = 2048;
pos_ARR = 1;
matrix_size = 64;

while matrix_size < max_matrix_size
    A = rand(matrix_size);
    [ A time_no  ] = BLAS3LU(A,64);
     [ A time_pp ] = BLAS3LUPP(A,64);
     
    time_NO_PP( pos_ARR ) = time_no;
    time_PP( pos_ARR ) = time_pp;
    
    matrix_dimension ( pos_ARR ) = matrix_size;
    matrix_size = matrix_size * 2;
    pos_ARR = pos_ARR +1;
end


plot(matrix_dimension,time_NO_PP,'s-','Color', color0);
hold on;
plot(matrix_dimension,time_PP,'d-','Color', color1);
hold on;

ylim([0 max(time_NO_PP*1.25)]);


l = legend('Time for solution without PP','Time for solution with PP' );


  




set(l,'FontSize',12);
ylabel('Bandwidth in Mbps');

xlabel('Time in sec.');
t = title({'Time for solution without PP','Time for solution with PP'},'interpreter','latex');

set(t,'FontSize',30);
set(gca,'fontsize',12);
