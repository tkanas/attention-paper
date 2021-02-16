%% Solve systems of equations for g vectors

N=length(covmatatt); % number of neurons
ncv=N*(N-1)/2; % number of cov values

eqind=zeros(ncv,4); % system of equations:
% columns represent [<neuron i>,<neuron j>,<unattended cov value>,<attended
% cov value>] as the equation CovA_ij=g_ig_jCovU_ij
ct=0;
for i=1:N
    for j=i+1:N
        ct=ct+1;
        eqind(ct,:)=[i,j,covmatunatt(i,j),covmatatt(i,j)];
    end
end

input0=rand(N,1); % Initialize solver
covanonf_of=@(inputs)covgivenpars_of(inputs,eqind,N); % Objective function to be minimized
opts=optimset('GradObj','on','Display','off');
[g,f] = fminunc(covanonf_of,input0,opts);
% g is the desired vector of g-values
