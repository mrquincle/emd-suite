#!/usr/bin/octave -qf

if (nargin < 3) 
	printf("Usage: display-file <cloud1> <cloud2> <match>\n")
	return
end

arg_list = argv ();

cloud1=arg_list{1};
cloud2=arg_list{2};
match=arg_list{3};

p1=load(cloud1);
p2=load(cloud2);
m=load(match);

% shift second cloud to origin
%M = size(p2,1)/2;
%p2a = p2(1:M,:);
%p2b = p2(M+1:M*2,:);
%p2a = p2a - mean(p2a);
%p2b = p2b - mean(p2b);
%p2 = [p2a; p2b];

%lim=5;
%p1=p1(1:lim,:);
%p2=p2(1:lim,:);

graphics_toolkit gnuplot

f = figure('visible','off')

m=m';

p=[p1; p2];

clf
%plot(p(:,1),p(:,2),'o');
plot(p1(:,1),p1(:,2),'ro');
axis('equal');
hold on
plot(p2(:,1),p2(:,2),'go');

[x1 x2] = meshgrid(p1(:,1),p2(:,1));
[y1 y2] = meshgrid(p1(:,2),p2(:,2));

% complete 
k=prod(size(x1));
%x1h = x1(:)(1:k);
%x2h = x2(:)(1:k);
%y1h = y1(:)(1:k);
%y2h = y2(:)(1:k);

x1h = x1(:);
x2h = x2(:);
y1h = y1(:);
y2h = y2(:);

xh=[x1h x2h];
yh=[y1h y2h];

% shift to origin
%xh = xh-mean(xh);
%yh = yh-mean(yh);

mflat=m(:);
mflat_sorted=sort(mflat,'descend');

% I can make a cool animation about this... I start with the largest streams and then gradually move to weights
% that are lower and lower and you'll see how first the ones that are furthest contribute most to the transportation
% plan and only later the ones that are closer to each other

% Hypothesis, epsilon-scaling with wormholes will first do the ones that are furthest apart (the corners)

m0=mflat_sorted(10);
m1=mflat_sorted(40);
m2=mflat_sorted(100);

idx=1:size(mflat,1);
idx=idx(mflat>=1)

for i=1:k
    if mflat(i)>=m0
	plot(xh(i,:), yh(i,:), '-b', 'linewidth', 2*log(mflat(i)*10+1e-7));
    elseif mflat(i)>m1
	plot(xh(i,:), yh(i,:), '-g', 'linewidth', 2*log(mflat(i)*10+1e-7));
    elseif mflat(i)>mflat_sorted(200)
	plot(xh(i,:), yh(i,:), '-c', 'linewidth', 2*log(mflat(i)*10+1e-7));
%   elseif mflat(i)>mflat_sorted(400)
%	plot(xh(i,:), yh(i,:), '-m', 'linewidth', 2*log(mflat(i)*10+1e-7));
  %  elseif mflat(i)>mflat_sorted(800)
		% do nothing
   % elseif mflat(i)>mflat_sorted(1000)
%	plot(xh(i,:), yh(i,:), '-m', 'linewidth', 2*log(mflat(i)*10));
    end
end
%plot([x1h x2h], [y1h y2h], 'o-r');

%line([x1(:) x2(:)], [y1(:) y2(:)])

print(f, "match.png", "-dpng")
