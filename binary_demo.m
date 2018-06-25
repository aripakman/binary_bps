% Author: Ari Pakman

% We compare two types of binary samplers:
% 1) Exact Hamiltonian Monte Carlo (HMC) for binary distributions
% 2) Binary Bouncy Particle Sampler (BPS)

% We use C++ implementations for Markov Random Fields (MRFs)
% In HMC we use the Gaussian augmentation and in BPS both Gaussian and
% Exponental augmentations

% The log probability of the MRF is defined as 
%
% log p(S|r,M) = -sum_i r(i)*S(i) - sum_{i<j} M(i,j)*S(i)*S(j)  + const  
%
% with S(i) = +1 or -1. 



%%
clear;

d = 10;    % dimension of binary vector 
frust = false;
c1 = .1;  % bias scale
c2 = .3;  % weights scale 
mrf = MRF(d,c1,c2, frust);    % create Markov Random Field object

M = mrf.M;
r = mrf.r;


%% Compile and run the HMC sampler 
mex  GCC='g++-4.9' ...
     COMPFLAGS='$COMPFLAGS -Wall -std=c++11'  ...
     src_cpp/hmc_binary.cpp ...
     src_cpp/MRF.cpp ...
     src_cpp/HMC_BinarySampler.cpp ...
     src_cpp/BinaryDistribution.cpp ...
 


% Warmup and discard 
L = 5000;
P = 50;
last_y = rand(1,d);
[XH,llH] = hmc_binary(M,r,L,P,last_y);


% Sample
L = 40000;
tic
[XH,llH] = hmc_binary(M,r,L,P,last_y);
cpu_time_hmc = toc;

mXH=mean(XH,1);


figure(99)
plot(mXH(1:1000))


% Signature of hmc_binary(...):

% Input
% M:        d x d symmetric matrix of real numbers
% r:        d x 1 vector of real numbers
% L:        number of samples
% P:        the travel time of the particle is (P+.5)*pi
% last_y :  d x 1 vector of real numbers. See explanation below
% seed:     optional input for random seed. If absent, the seed is chosen from the C rand() function

% Output 
% XH:   d X L matrix of samples 
% llH:  L-dim vector, each entry is the log-likelihood of each sample


%% Role of last_y

% The Markov chain is defined over the continuos variables y, such that
% S = sign(y). The sampler requires initial values for y in the
% variable last_y. When the function returns, last_y will contain the last 
% values of y (last_y is passed by reference). 
% This is useful when the binary variables are part of a Gibbs sampling scheme.



%% Compile and run the BPS sampler 
mex  GCC='g++-4.9'...
     COMPFLAGS='$COMPFLAGS -Wall -std=c++11'  ...
     src_cpp/bouncy_binary.cpp ...
     src_cpp/MRF.cpp ...
     src_cpp/BouncyBinarySamplerExponential.cpp ...
     src_cpp/BouncyBinarySamplerGaussian.cpp ...
     src_cpp/BinaryDistribution.cpp
 

 
 
max_time = 1e6;
last_y = rand(1,d)-.5;


% Exponential Augmentation
tic
[XE, timesE, llE] = bouncy_binary(0, M,r, max_time ,last_y);
cpu_time_bouncyE = toc


% Gaussian Augmentation
tic
[XG, timesG, llG ] = bouncy_binary(1, M,r, max_time ,last_y);
cpu_time_bouncyG = toc





% Signature of hmc_binary(...):

% Input
% Augmentaton type:     0 = Exponential, 1 = Gaussian
% M:                    d x d symmetric matrix of real numbers
% r:                    d x 1 vector of real numbers
% max_time:             maximum travel time of the particle
% last_y :              d x 1 vector of real numbers. See explanation below
% seed:                 optional input for random seed. If absent, the seed is chosen from the C rand() function

% Output 
% X:     d X LB matrix of samples 
% times: LB-dim vector, each entry is the time spent in each on the LB sample 
% ll:    LB-dim vector, each entry is the log-likelihood of each sample

% the number of different samples returned, LB, 
% depends on the max_times alloted to the travelling particle



%%



% histogram of travel times for the bouncy samplers
figure(99)
hist(timesE, 160)


% plot loglikelihoods of the HMC samples
figure(12)
clf()
hold on
plot(llH)
grid()


%% transform the bouncy samples to discrete samples in order to compare the ACF with HMC
% we take into account the CPU time

time_per_sample_hmc = cpu_time_hmc/L;

% exponential
LbE = floor(cpu_time_bouncyE/time_per_sample_hmc);  % number_of_bouncy_samples  

dt = sum(timesE)/LbE;      
llsE = zeros(LbE,1);
mXE = zeros(LbE,1);

i = 1;
j = 1;
tt =0;

while j < LbE +1
     
     tt = tt +timesE(i);
     if tt > dt 
         llsE(j) = llE(i);
         mXE(j) = mean(XE(:,i));
         j = j+ 1;   
         tt = tt-dt;
     end
    i = i +1; 
    if i > size(XE,2)
        break
    end
    
end

% Gaussian
Lb = floor(cpu_time_bouncyG/time_per_sample_hmc);  % number_of_bouncy_samples  

dt = sum(timesG)/Lb;      
llsG = zeros(Lb,1);
mXG = zeros(Lb,1);
i = 1;
j = 1;
tt =0;

while j < Lb +1
     
     tt = tt +timesG(i);
     if tt > dt 
         llsG(j) = llG(i);
         mXG(j) = mean(XG(:,i));
         
         j = j+ 1;   
         
         tt = tt-dt;
     end
    i = i +1; 
    if i > size(XG,2)
        break
    end
    
end





%% plot the autocorrelation functions 

fig=figure(66);
clf
hold on
plot(acf(llsE,10))
plot(acf(llsG,10))
plot(acf(llH,10))
grid
title('ACF')
xlabel('lag')
axis tight
legend('BPS Exponential', 'BPS Gaussian', 'HMC')
ylim([-0.1,1])
box on
set(gca,'XTickLabel', 0:10)





%%  plot magnetization and log probabilities for first 1000 samples

figure(10)
clf();
hold on
m = mean(XG,1);
plot(m(1:3000))

mh = mean(XH,1);
plot(mh(1:3000))
legend('BPS Gaussian', 'HMC')






