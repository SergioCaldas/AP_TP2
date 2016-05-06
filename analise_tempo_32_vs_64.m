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
    
    [ A_3pi  time_3pp error3PIVOT ] = BLAS3LUPP(A,64);
    [ A_2pi  L_2pi U_2pi P_2pi time_2pp error2PIVOT ] = BLAS2LUPP(A);

    time_2NO_PP_64( pos_ARR ) = time_2no;
    time_2PP_64( pos_ARR ) = time_2pp;
    
    time_3NO_PP_64( pos_ARR ) = time_3no;
    time_3PP_64( pos_ARR ) = time_3pp;
    
    matrix_dimension ( pos_ARR ) = matrix_size;
    matrix_size = matrix_size * 2;
    pos_ARR = pos_ARR +1;
end

matrix_size = 64;
pos_ARR = 1 ;

while matrix_size <= max_matrix_size
    A = rand(matrix_size);
    [ A_3no   L_3no U_3no time_3no error3NOPIVOT  ] = BLAS3LU(A,32);
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
name = {'compute 431';'compute 641';'compute 652';'compute 662'};
%x = [ 64 128 256 512 1024 ]; 

x = [1:5]; 

bpcombined = [time_2NO_PP(:), time_2NO_PP_64(:), time_3NO_PP(:) , time_3NO_PP_64(:) , time_2PP(:), time_2PP_64(:), time_3PP(:) , time_3PP_64(:)];
hb = bar(x, bpcombined, 'grouped')

bg = [1 1 1; 0 0 0]
cores = distinguishable_colors(27,bg)
figure(1)
set(hb(1), 'FaceColor',cores(1,:))
set(hb(2), 'FaceColor',cores(2,:))
set(hb(3), 'FaceColor',cores(3,:))
set(hb(4), 'FaceColor',cores(4,:))
set(hb(5), 'FaceColor',cores(5,:))
set(hb(6), 'FaceColor',cores(6,:))
set(hb(7), 'FaceColor',cores(7,:))
set(hb(8), 'FaceColor',cores(8,:))

%ylim([0 max(time_3PP*2)]);
%xlim([0 max(matrix_dimension)]);


l = legend('Time for solution without PP (block size 32) -- BLAS2LU', 'Time for solution without PP (block size 64) -- BLAS2LU',  'Time for solution without PP (block size 32) -- BLAS3LU', 'Time for solution without PP (block size 64) -- BLAS3LU',  'Time for solution with PP (block size 32) -- BLAS2LUPP', 'Time for solution with PP (block size 64) -- BLAS2LUPP',  'Time for solution with PP (block size 32) -- BLAS3LUPP', 'Time for solution with PP (block size 64) -- BLAS3LUPP' );
set(l,'Location','southoutside')

set(l,'FontSize',12);
ylabel('Time in sec.');
%set(gca,'Xtick',matrix_dimension);

xlabel('Matrix Dimension');
t = title({'BLAS2LU, BLAS2LUPP, BLAS3LU, BLAS3LUPP time for solution analysis ','based on matrix dimension, for block size 32 and 64'},'interpreter','latex');

set(t,'FontSize',30);
set(gca,'fontsize',12);

set(gca,'fontsize',12);
name = {'64';'128';'256';'512'; '1024'};

set(gca,'xticklabel',name);
