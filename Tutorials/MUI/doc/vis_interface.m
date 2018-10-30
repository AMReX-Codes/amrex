clear; clc; close all;

graphics_toolkit fltk;
%graphics_toolkit gnuplot;

N3D = [4 4 4];
nstep3D = [2 2 2];

N2D = N3D(1:2);
nstep2D = [1 1];

sizes3D = fix(N3D./nstep3D)+1;
x3d = (0:nstep3D(1):N3D(1))'; y3d = (0:nstep3D(2):N3D(2))'; z3d = (0:nstep3D(3):N3D(3))';
[X3D,Y3D,Z3D] = meshgrid(x3d,y3d,z3d);

% Rectangular
sizes2D = fix(N2D./nstep2D)+1;
x2d = (0:nstep2D(1):N2D(1))'; y2d = (0:nstep2D(2):N2D(2))';
[X2D,Y2D] = meshgrid(x2d,y2d);

%% Figure properties

grid_fig = figure;
hold on

view(35,35)
axis equal
axisfontsize = 14;
xlabel('x','FontSize',axisfontsize)
ylabel('y','FontSize',axisfontsize)
zlabel('z','FontSize',axisfontsize)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])


%% Plot 3D gridlines

for k = 1:sizes3D(3)
    plot3(X3D(:,:,k),Y3D(:,:,k),Z3D(:,:,k),'--k','Linewidth',1);
end
for k = 1:sizes3D(3)
    plot3(X3D(:,:,k)',Y3D(:,:,k)',Z3D(:,:,k)','--k','Linewidth',1);
end
for j = 1:sizes3D(2)
    plot3(reshape(X3D(j,:,:),sizes3D([1 3]))',reshape(Y3D(j,:,:),sizes3D([1 3]))', ...
    reshape(Z3D(j,:,:),sizes3D([1 3]))','--k','Linewidth',1);
end

%% Plot 2D gridlines

% Rectangular
plot(X2D,Y2D,'r-','Linewidth',2);
plot(X2D',Y2D','r-','Linewidth',2);

%% Draw 3D faces (for MATLAB)

% rng(777);
% COL = rand(prod(N3D),3);
% 
% i_p = ones(4,1); j_p = i_p; k_p = i_p;
% f = [1 2 3 4];
% 
% index = 1;
% for k = 1:sizes3D(3)-1
% for j = 1:sizes3D(2)-1
% for i = 1:sizes3D(1)-1
%     for facedir = 1:3
%         for hilo = 0:1
%             i_p(:) = i; j_p(:) = j; k_p(:) = k;
%             switch facedir
%                case 1
%                   i_p = i_p + hilo;
%                   j_p = j_p + [0 0 1 1]'; 
%                   k_p = k_p + [0 1 1 0]';
%                case 2
%                   j_p = j_p + hilo;
%                   k_p = k_p + [0 0 1 1]'; 
%                   i_p = i_p + [0 1 1 0]';
%                case 3
%                   k_p = k_p + hilo;
%                   i_p = i_p + [0 0 1 1]'; 
%                   j_p = j_p + [0 1 1 0]';
%             end
%             v = [x3d(i_p) y3d(j_p) z3d(k_p)];
%             col = COL(index,:);
%             patch('Faces',f,'Vertices',v,'FaceColor',col, ... 
%                 'EdgeColor','none','FaceAlpha',.2);
%         end
%     end
%     index = index+1;
% end
% end
% end

%% Save figure
%saveas(grid_fig,'grid.png');