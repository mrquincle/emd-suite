#!/usr/bin/octave -qf

if (nargin < 4) 
	printf("Usage: display-file <cloud1> <offset1> <cloud2> <offset2>\n")
	return
end

arg_list = argv ();

cloud1=arg_list{1};
offset1=arg_list{2};
cloud2=arg_list{3};
offset2=arg_list{4};

orig1=load(cloud1);
offset1=load(offset1);

pnts = orig1 - offset1;

graphics_toolkit gnuplot

fig=figure('visible','off')

plot(pnts(:,1),pnts(:,2),'go');
axis('equal');
hold on

orig2=load(cloud2);
offset2=load(offset2);

n=size(orig2,1);
k=n/2;
orig2(k:n,:);
offset2(k:n,:);

pnts = orig2 - offset2;

plot(pnts(:,1),pnts(:,2),'ro');

print(fig, "result.png", "-dpng")
