clear all;
clf;
set(gca,'FontSize',15);
set(0,'DefaultAxesFontSize', 22)

%------------------ Scaling;

load runtime.dat
x = runtime(:,1);
y = runtime(:,2);

fac = y(1);
y = fac./y;

x1 = [64 512 1024 2048 5096 100000];
y1 = [1 1 1 1 1 1];

h = semilogx(x,y,'bo-',x1,y1,'k:');
set(h,'LineWidth',2);
xlabel('Number of cores','FontSize',22);
ylabel('Efficiency','FontSize',22);
title('Weak scaling');
l = legend('WarpX','Ideal','Location','SouthEast');
l.FontSize = 18;
ax = axis;
axis([50 10000 0 1.2]);

print -dpng wscaling.png

pause;

%------------------ Profiling;

load profile.dat
ncore = profile(:,1);
arrayCopy = profile(:,2);
particleCopy = profile(:,3);
picsarCurrentDep = profile(:,4);
picsarPartPush = profile(:,5);
particleEvolve = profile(:,6);
redist = profile(:,7);
redistMPI = profile(:,8);
picsarFieldGather = profile(:,9);

h = semilogx(ncore,particleCopy,'bo-.', ...
	     ncore,arrayCopy,'ro-.', ...
	     ncore,particleEvolve,'ko-.', ...
	     ncore,redist,'go-.', ...
	     ncore,redistMPI,'mo-.', ...
	     ncore,picsarCurrentDep,'b*-', ...
	     ncore,picsarPartPush,'r*-', ...
	     ncore,picsarFieldGather,'k*-');
set(h,'LineWidth',2);
xlabel('Number of cores','FontSize',22);
ylabel('Maximum percentage','FontSize',22);
title('Loads of WarpX code components');
l = legend('Particle copy', ...
	   'Array copy', ...
	   'Particle evolve', ...
	   'Particle redistribute', ...
	   'Redistribute MPI', ...
	   'Current deposit', ...
	   'Particle push', ...
	   'Field gather', ...
	   'Location','NorthEast');
l.FontSize = 13;
ax = axis;
axis([50 10000 0 50]);

print -dpng profile.png


