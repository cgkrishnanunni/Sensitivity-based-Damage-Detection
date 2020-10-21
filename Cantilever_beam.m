
% Programmed by Krishnanunni CG, NIT CALICUT-2017%
% Refer Sensitivity-based damage detection algorithm for structures 
% using vibration data. J Civil Struct Health Monit 9, 137–151 (2019). 
% https://doi.org/10.1007/s13349-018-0317-0

function [bestnest,fmin]=cantilever()

clc;
clear;
L=0.4953;  %input('Length of beam');
b=0.0254;  %input('Breadth of beam');
d=0.00635; %input('Depth of beam');
nElems=20; %input('Enter number of Elements');
nNodes=nElems+1;
Ar=b*d;
E=7.1e10;  %input('Modulus of Elasticity');
p=2210;    %input('Density ');
I=(b*d*d*d)/12;
gStiff=zeros((nElems*3),(nElems*3));
gMass=zeros((nElems*3),(nElems*3));
dgStiff=zeros((nElems*3),(nElems*3));
Ncrd=zeros(2,nNodes);
elemConn=zeros(2,nElems);
%Nodal Coordinates
for i=1:nNodes
   Ncrd(1,i)=(L/nElems)*(i-1);
   Ncrd(2,i)=0;   
end


%Element connectivity

for i=1:nElems
   elemConn(1,i)=i;
   elemConn(2,i)=i+1;    
end


%Rest Node Condition
Rnod=zeros(3,nNodes);
Restnodes=1; %input('enter the number of rest nodes=');
Rnod(1,1)=1;
Rnod(2,1)=1;
Rnod(3,1)=1;

for j=2:nNodes
    for i=1:3
        Rnod(i,j)=0;
    end
end
 
neq=0;

for j=1:nNodes
    for i=1:3
        if (Rnod(i,j)==0)
            neq = neq+1;
				Rnod(i,j)=neq;
       		continue;	
        else
				Rnod(i,j)=0;        
        end
    end
end  



%Element Stiffness
Le=L/nElems;%element length
A=(Ar*E)/Le;
B=(12*E*I)/(Le^3);
C=(6*E*I)/(Le^2); 
D=(4*E*I)/Le;
F=(2*E*I)/Le;
eStiff=[A 0 0 -A 0 0;0 B C 0 -B C;0 C D 0 -C F;-A 0 0 A 0 0;0 -B -C 0 B -C;0 C F 0 -C D];


%Element Mass matrix
m=(p*Ar*Le)/420;
G=22*Le;
H=13*Le;
J=4*Le*Le;
Q=3*Le*Le;
eMass=m*[140 0 0 70 0 0;0 156 G 0 54 -H;0 G J 0 H -Q;70 0 0 140 0 0;0 54 H 0 156 -G;0 -H -Q 0 -G J];


%Assembly of element Matrices
kk=zeros(1,6);
xl=zeros(2,1);
yl=zeros(2,1);

for n=1:nElems
    for i=1:2        
        xl(i)=Ncrd(1,elemConn(i,n));
        yl(i)=Ncrd(2,elemConn(i,n));
    end
    dof=0;   
    for i=1:2
        for j=1:3
            dof=dof+1;
            kk(dof)=Rnod(j,elemConn(i,n));
        end
    end
   % display(kk)    
   for i=1:6
           if(kk(i)==0)
           continue;
           end
           k=kk(i);
          
         for j=1:6
           if(kk(j)==0)
           continue;
           end
           gStiff(k,kk(j))=gStiff(k,kk(j))+eStiff(i,j);
           gMass(k,kk(j))=gMass(k,kk(j))+eMass(i,j);
         end
   end
   

end     


 W=inv(gMass)*gStiff; 
 w=zeros(nElems*3);
 w=eig(W); % Eigen values of the beam
 

 [T,X]=eig(W);
 display(T)


if nargin<1,
% Number of nests (or different solutions)
n=25;
end

% Discovery rate of alien eggs/solutions
pa=0.35;

%% Change this if you want to get better results
% Tolerance
Tol=0.01;
%% Simple bounds of the search domain
% Lower bounds
nd=nElems; 
Lb=0*ones(1,nd); 
% Upper bounds
Ub=1*ones(1,nd);

% Random initial solutions
for i=1:n,
nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end

% Get the current best
fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w);

N_iter=0;
%% Starting iterations
while (N_iter<65000),
    

    % Generate new solutions (but keep the current best)
     new_nest=get_cuckoos(nest,bestnest,Lb,Ub);   
     [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w);
    % Update the counter
      N_iter=N_iter+n; 
    % Discovery and randomization
      new_nest=empty_nests(nest,Lb,Ub,pa) ;
    
    % Evaluate this set of solutions
      [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w);
    % Update the counter again
      N_iter=N_iter+n;
    % Find the best objective so far  
    if fnew<fmin,
        fmin=fnew;
        bestnest=best;
    end
    
end %% End of iterations


%% Post-optimization processing
%% Display all the nests
disp(strcat('Total number of iterations=',num2str(N_iter)));
fmin
bestnest

%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best,Lb,Ub)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
    s=nest(j,:);
    % This is a simple way of implementing Levy flights
    % For standard random walks, use step=1;
    %% Levy flights by Mantegna's algorithm
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
  
    % In the next equation, the difference factor (s-best) means that 
    % when the solution is the best solution, it remains unchanged.     
    stepsize=0.01*step.*(s-best);
    % Here the factor 0.01 comes from the fact that L/100 should the typical
    % step size of walks/flights where L is the typical lenghtscale; 
    % otherwise, Levy flights may become too aggresive/efficient, 
    % which makes new solutions (even) jump out side of the design domain 
    % (and thus wasting evaluations).
    % Now the actual random walks or flights
    s=s+stepsize.*randn(size(s));
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);
end


%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w)
% Evaluating all new solutions
for j=1:size(nest,1),
    fnew=fobj(newnest(j,:),nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w);
    if fnew<=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end
% Find the current best
[fmin,K]=min(fitness) ;
best=nest(K,:);
%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
  new_nest(j,:)=simplebounds(s,Lb,Ub);  
end

% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;

%% You can replace the following by your own functions
% A d-dimensional objective function
function z=fobj(u,nElems,elemConn,Ncrd,Rnod,gStiff,eStiff,gMass,eMass,dgStiff,T,w)



%Assembly of Damaged element Matrices
kk=zeros(1,6);
xl=zeros(2,1);
yl=zeros(2,1);
dgStiff=zeros((nElems*3),(nElems*3));
for n=1:nElems
    for i=1:2        
        xl(i)=Ncrd(1,elemConn(i,n));
        yl(i)=Ncrd(2,elemConn(i,n));
    end
    dof=0;   
    for i=1:2
        for j=1:3
            dof=dof+1;
            kk(dof)=Rnod(j,elemConn(i,n));
        end
    end
   % display(kk)    
   for i=1:6
           if(kk(i)==0)
           continue;
           end
           k=kk(i);
          
         for j=1:6
           if(kk(j)==0)
           continue;
           end
           dgStiff(k,kk(j))=dgStiff(k,kk(j))+(u(n)*eStiff(i,j));
           
         end
   end
   

end     
dgStiff=-dgStiff;


% EXPERIMENTAL DATA FOR SINGLE CRACK (EXAMPLE 1)

exp1=[19.53;122.05;339.26;661.73;1085.22;1594.59]; % Undamaged
exp2=[19;115.85;332.36;646.91;1037.46;1591.36];    % Damaged

for i=1:6
    exp1_omega(i)=4*pi*pi*exp1(i)*exp1(i);
    exp2_omega(i)=4*pi*pi*exp2(i)*exp2(i);
end

for i=1:6
vd(i)=(exp2_omega(i)-exp1_omega(i))/exp1_omega(i);
vd_new(i)=(exp2(i)-exp1(i))/exp1(i);
end
f=0;


for i=1:6
    
    dlamda(i,1)=(transpose(T(:,(nElems*3)-i+1))*dgStiff*T(:,(nElems*3)-i+1))/(transpose(T(:,(nElems*3)-i+1))*gMass*T(:,(nElems*3)-i+1));
    
    f=f+((dlamda(i,1)/w((nElems*3)-i+1))-vd(i))^2;
    
end

% ALTERNATE OBJECTIVE FUNCTION USED BY PREVIOUS RESEARCHERS
% COMPUTATIONALLY EXPENSIVE OBJECTIVE FUNCTION DUE TO EIGEN VALUE
% COMPUTATION AND REQUIRES MORE ITERATIONS

 %gStiff_damaged=gStiff+dgStiff;
 %W_damaged=inv(gMass)*gStiff_damaged; 
 %w_damaged=zeros(nElems*3);
 %w_damaged=sqrt(eig(W_damaged))/(2*pi); % Eigen values of the beam
 %w=sqrt(w)/(2*pi);
 
 %for i=1:6
    
    
    %f=f+(((w_damaged((nElems*3)-i+1)-w((nElems*3)-i+1))/w((nElems*3)-i+1))-vd_new(i))^2;
    
%end
 
 
f=sqrt(f) % VALUE OF Objective function
z=f;
