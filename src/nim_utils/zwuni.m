% computes the N+1 uniform nodes z and weights on [-1,1]
function [z,w] = zwuni(N)

n = N+1;

z=(0:N)'; z = -1 + 2*z./N;

w=ones(n,1); w=2*w./N;
w(1)=.5*w(1);
w(n)=.5*w(n);
end