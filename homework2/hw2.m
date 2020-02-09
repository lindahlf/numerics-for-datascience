%% Part 2a 
% Show all steps of k-means explicitly for a simple example
clear all; close all; clc

X = [1 1 1; 1.5 2 0; -0.5 0 -3; 0 -0.5 -1]; % datapoints
R = [3 3 3; -1 -1 0]; % starting centers 

% First iteration
Dist = dist(X,R')'; % Computing distance between X and R 

IA = [0 1 0 0]'; IB = [1 0 1 1]';

V = [normalize(IA,'norm', 1) normalize(IB,'norm', 1)]
R = V'*X % new centers 

% Second interation
Dist = dist(X,R')'

IA = [1 1 0 0]'; IB = [0 0 1 1]';

V = [normalize(IA,'norm', 1) normalize(IB,'norm', 1)]
R = V'*X % new centers 

% Third iteration
Dist = dist(X,R')'

%% Part 2b
% Show all steps of k-means explicitly for another simple example
clear all; close all; clc

X = [1 1 1; 1.5 2 0; -0.5 0 -3; 0 -0.5 -1]; % datapoints
R = [-0.5 0 -3; -1 -1 1]; % starting centers 

% First iteration
Dist = dist(X,R')' % Computing distance between X and R 

IA = [0 0 1 1]'; IB = [1 1 0 0]';

V = [normalize(IA,'norm', 1) normalize(IB,'norm', 1)]
R = V'*X % new centers 

% Second interation
Dist = dist(X,R')'

%% Part 3 
% Building a similarity graph
clear all; close all; clc

ee = 2.5; % epsilon 

n0 = 5; p = 3;
randn('seed',0);
c = 2;
A1 = randn(n0,p)+c*[1,0,0];
A2 = randn(n0,p)+c*[0,1,0];
A3 = randn(n0,p)+c*[0,0,1];
A4 = [0,0,0];
A = [A1;A2;A3;A4]; % Randomly generated matrix

Dist = dist(A,A'); % distance matrix 

% Constructing the weight matrix of an epsilon-neighborhood graph from the distance matrix.
W = ones(size(Dist));
W(Dist > ee) = 0; % Set too far away to
W(Dist == 0) = 0;  % Avoid an edge to itself


% Constructing the degree matrix
D = zeros(size(A,1));

for i=1:size(A,1)
    D(i,i)=sum(W(i,:));
end

% Plotting the graph from the weight matrix

n = size(W,1);
z = exp(((1:n)-1)*2i*pi/n);
clf;

plot(z,'o','MarkerFaceColor','k'); % Plot all nodes
hold on;

for i=1:size(W,1)
   for j=1:size(W,1)
       if (W(i,j)>0)
            plot([z(i),z(j)],'k');   % Plot edges
       end
   end
end

% b) 
close all
L = D - W

%% Part 3c-d
clear all; close all; clc

ee = 2.5; % epsilon 

n0 = 5; p = 3;
randn('seed',0);
c = 3;
A1 = randn(n0,p)+c*[1,0,0];
A2 = randn(n0,p)+c*[0,1,0];
A3 = randn(n0,p)+c*[0,0,1];
A4 = [0,0,0];
A = [A1;A2;A3;A4]; % Randomly generated matrix

Dist = dist(A,A'); % distance matrix 

% Constructing the weight matrix of an epsilon-neighborhood graph from the distance matrix.
W = ones(size(Dist));
W(Dist > ee) = 0; % Set too far away to be connected
W(Dist == 0) = 0;  % Avoid an edge to itself


% Constructing the degree matrix
D = zeros(size(A,1));

for i=1:size(A,1)
    D(i,i)=sum(W(i,:));
end

% Plotting the graph from the weight matrix

n = size(W,1);
z = exp(((1:n)-1)*2i*pi/n);
clf;

plot(z,'o','MarkerFaceColor','k'); % Plot all nodes
hold on;

for i=1:size(W,1)
   for j=1:size(W,1)
       if (W(i,j)>0)
            plot([z(i),z(j)],'k');   % Plot edges
       end
   end
end

L = D - W;
eig(L)

%% Part 5
clear all; close all; clc
load('bengali_cleanup.mat');

A = imread('bengali_map.png');
figure(1); clf;
imshow(A); hold on;
jv = [102,280,10];
%plot(y_coords(jv), x_coords(jv),'r*') % Note: x and y reversed since images have swapped x and y axis.
plot(y_coords, x_coords,'r*')
%title("Measurement stations 102, 280 and 10")
title("All measurement stations")
figure(2);
plot(tv, timeseries(jv(1),:), '-'); hold on
plot(tv, timeseries(jv(2),:), '--')
plot(tv, timeseries(jv(3),:), '-.')


Dist = pdist2(timeseries, timeseries);
Dist(102,280)
Dist(102,10)
Dist(280,10)

%% Part 5c 
% k-nearest neighbours 
clear all; close all; clc
load('bengali_cleanup.mat');

k = 3;

Dist = pdist2(timeseries, timeseries);
W = zeros(size(Dist));


% Finds k closest neighbours 
[D,I] = pdist2(timeseries,timeseries,'euclidean','Smallest',k+1); 
D = D(2:k+1,:);
I = I(2:k+1,:);

% kNN (version: or)
for i = 1:937
    W(i,I(:,i)) = 1;
    W(I(:,i),i) = 1;
end

% Loop to find weight matrix, kNN (version: and) 
% for i = 1:length(x_coords)
%     temp = I == i;
%     [m,n] = find(temp);
%     for elem = n
%         if any(ismember(I(:,i),n)) == 1
%             W(i,n) = 1;
%             W(n,i) = 1;
%         end
%     end
% end

%clf; plot(graph(W));




