%The computational region is a square of size xlength*ylength, i.e. Dim(1)*Dim(2), with pml of
%'thickness' enclosed both in x and y.

%Any function (source, epsilon, solution) is defined on the M*N grid points

%The ordering of the grid goes like [1 2 ... N; N+1 N+2 ... 2N; ...] 

%This is an eigenproblem.
%Eigen_Maxwell is the eps^(-1)curl mu^(-1) curl operator.
%The user can obtain eigenmodes by eigensolvers in Matlab

%pml size and strength should be chosen carefully to unveil QNMs and the
%continuum in a distinctive manner.


addpath(genpath('TE'));
format long;
%% setups
h = 0.01; %size of a single grid 
omega = 2*pi/1; %%resonant frequecy 
design_region = [3.4;3.4]; %% [x-length;y-length]
pml_thickness = [0.5;0.5]; %% thickness of pml [x;y]
dim = design_region + 2.*pml_thickness;
N = round(dim(1) / h) + 1; %num of x dim grid points
M = round(dim(2) / h) + 1; %num of y dim grid points

beta = 20; %%pml strength
BC = {{'pml', [pml_thickness(1), beta]}, {'pml', [pml_thickness(2), beta]}}; %%boundry condition {x,y}

% epsilon = ones(M,N); % one can define the epsilon here, material should be away from PML
eps = readmatrix('test_TE.txt'); 
imagesc(eps);
[szy,szx] = size(eps);
epsilon = ones(M,N);
epsilon(round((M-szy)/2)+1:round((M-szy)/2)+szy,round((M-szy)/2)+1:round((M-szy)/2)+szy) = eps;

%% find modes and plot
num_mode = 10; %% we will find # of modes with frequncy closest to the pre-defined omega
[rigeig,rigeigval] = ModeSolverTE(h,dim,BC,epsilon,num_mode,omega);

figure;
subplot(1,2,1);imagesc(real(reshape(rigeig(:,1),N,M).'));
subplot(1,2,2);imagesc(imag(reshape(rigeig(:,1),N,M).'));
omega_mode = sqrt(rigeigval(1,1));
display(['omega is ' num2str(omega_mode)]);
