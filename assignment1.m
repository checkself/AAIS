tsx=logRets(:,1);
cac=logRets(:,2);
dax=logRets(:,3);
eurostox=logRets(:,4);
nikkei=logRets(:,5);
ftse=logRets(:,6);
sp=logRets(:,7);
ibov=logRets(:,8);
N=size(logRets,2);

%% Question1

means=mean(logRets);
stds=std(logRets);
skew=skewness(logRets);
kurt=kurtosis(logRets);
corrMat=corr(logRets);
covMat=cov(logRets);


%figure;
for i=1:N;
    % using J-B test because it tests for deviations in the 3rd and 4th
    % moments from the normal
   [jb.h(i),jb.p(i),jb.jbstat(i),jb.critval(i)]=jbtest(logRets(:,1));
   %subplot(2,4,i);
   %normplot(logRets(:,i));
   %title(indecies{i});
end


%% Question2

H = 2.0 * covMat;
f = zeros(N,1);
% 1st Constraint: weights sum to one
Aeq = ones(1,N);
beq = 1.0;
% 2nd Constraint: short sales
A = zeros(N);
b = zeros(N,1);
[weights1, sigma1] = quadprog(H, f, A, b, Aeq,beq, [], [], [], []);


% 2nd Constraint: no short sales
A = -eye(N);
b = zeros(N,1);
[weights2, sigma2] = quadprog(H, f, A, b, Aeq,beq, [], [], [], []);

%% Question 3
d=dates(1:end-13);
windowSize=5;
inRange=zeros(length(d),length(breakDates)-windowSize);
for i=1:length(breakDates)-windowSize
    inRange(:,i)=(d>breakDates(i))&(d<breakDates(i+windowSize));
    rows=find(inRange(:,i));
    rollingLogRets{i}=logRets(rows,:);
end

