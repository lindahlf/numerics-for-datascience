%% Part 2a
clear all; close all; clc

b = [1 1 2 2]';

w = @(N) exp(-2*pi*i/N);

bhat = fft(b)

% bajs på den här frågan

%% Part 2b: implement the FFT 
clear all; close all; clc

b = tan(2*(1:16)');

fftx(b)

%% Part 2c: compare own FFT with naive-matrix version
clear all; close all; clc

nv = 2.^(1:13);

P=50; % number of samples

% Note: It is important to preallocate timing matrices, and not dynamically expand them
T1=zeros(length(nv),P);  
T2=zeros(length(nv),P);

for p=1:P
    disp(p)
    for k=1:length(nv)
        n=nv(k);
        b=randn(n,1);

        % fftx.
        tic
        z=fftx(b);
        T1(k,p)=toc;

        % DFT 
        tic
        z=dftmtx(n)*b;
        T2(k,p)=toc;
    end
end

t1=mean(T1,2); % Plot the mean of the runs
t2=mean(T2,2);

loglog(nv,t1,nv,t2);
legend('FFT','naive');
xlabel('n');
ylabel('CPU time');

%% Part 2d: improved fft for symmetric vectors 
clear all; close all; clc 


nv = 2.^(1:13);

P=50; % number of samples

% Note: It is important to preallocate timing matrices, and not dynamically expand them
T1=zeros(length(nv),P);  
T2=zeros(length(nv),P);
T3=zeros(length(nv),P);


for p=1:P
    disp(p)
    for k=1:length(nv)
        n=nv(k);
        b1 = randn(n/2,1);
        b = reshape(repmat(b1(:).',2,1),1,[])';

        % fftx.
        tic
        z=fftx(b);
        T1(k,p)=toc;

        % DFT 
        tic
        z=dftmtx(n)*b;
        T2(k,p)=toc;
        
        % FFTsym
        tic 
        z = fftsym(b,1);
        T3(k,p)=toc;
    end
end

t1=mean(T1,2); % Plot the mean of the runs
t2=mean(T2,2);
t3=mean(T3,2);


loglog(nv,t1,nv,t2,nv,t3);
legend('FFT','naive','Symmetric');
xlabel('n');
ylabel('CPU time');

%% Part 3a: audio
clear all; close all; clc

[Y, FS]=audioread("hw3_terrible_sound_with_hidden_message.ogg.");

yhat = fft(Y);

%plot((1:1:length(yhat)), abs(yhat))

I = find(abs(yhat)>18);
yhat(I) = 0;

goodsound = ifft(yhat);
% Numerics is fun 

%% Part 3b: large problems and DFT
clear all; close all; clc


[Y, FS]=audioread("hw3_terrible_sound_with_hidden_message.ogg.");

F = dftmtx(length(Y));

yhat = F*Y;

%% Part 4
clear all; close all; clc
% Compute a semiseparable matrix of size n times n
n=500;
randn('seed',0);
v=randn(n,1);
D=diag(randn(n,1));
A=v*v'+D;


    









