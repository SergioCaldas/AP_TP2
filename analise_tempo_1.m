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


max_matrix_size = 1024;
pos_ARR = 1;
matrix_size = 64;

while matrix_size <= max_matrix_size
    A = rand(matrix_size);
    [ A_3no   L_3no U_3no time_3no error3NOPIVOT  ] = BLAS3LU(A,64);
    [ A_2no   L_2no U_2no time_2no error2NOPIVOT ] = BLAS2LU(A);
    
    [ A_3pi  time_3pp error3PIVOT ] = BLAS3LUPP(A,32);
    [ A_2pi  L_2pi U_2pi P_2pi time_2pp error2PIVOT ] = BLAS2LUPP(A);

    time_2NO_PP( pos_ARR ) = time_2no;
    time_2PP( pos_ARR ) = time_2pp;
    
    time_3NO_PP( pos_ARR ) = time_3no;
    time_3PP( pos_ARR ) = time_3pp;
    
    matrix_dimension ( pos_ARR ) = matrix_size;
    matrix_size = matrix_size * 2;
    pos_ARR = pos_ARR +1;
end

FigHandle = figure;
set(FigHandle, 'Position', [0, 0, 640, 480]);
max_matrix_size = 2048;
pos_ARR = 1;
matrix_size = 64;

plot(matrix_dimension,time_2NO_PP,'s-','Color', color0);
hold on;
plot(matrix_dimension,time_3NO_PP,'s-','Color', color1);
hold on;
plot(matrix_dimension,time_2PP,'d-','Color', color2);
hold on;
plot(matrix_dimension,time_3PP,'d-','Color', color7);
hold on;


%ylim([0 max(time_3PP*2)]);
xlim([0 max(matrix_dimension)]);


l = legend('Time for solution without PP -- BLAS2LU', 'Time for solution without PP -- BLAS3LU', 'Time for solution with PP -- BLAS2LUPP', 'Time for solution with PP -- BLAS3LUPP' );

set(l,'FontSize',12);
ylabel('Time in sec.');
set(gca,'Xtick',matrix_dimension);

xlabel('Matrix Dimension');
t = title({'BLAS2LU, BLAS2LUPP, BLAS3LU, BLAS3LUPP time for solution analysis ','based on matrix dimension, for block size 64'},'interpreter','latex');

set(t,'FontSize',30);
set(gca,'fontsize',12);


