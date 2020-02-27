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
















