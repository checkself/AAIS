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
A_short = zeros(N);
b_short = zeros(N,1);
[weightsTotalShort, sigmaTotalShort] = quadprog(H, f, A_short, b_short, Aeq,beq, [], [], [], []);


% 2nd Constraint: no short sales
A_noShort = -eye(N);
b_noShort = zeros(N,1);
[weightsTotalNoShort, sigmaTotalNoShort] = quadprog(H, f, A_noShort, b_noShort, Aeq,beq, [], [], [], []);

%% Question 3
d=dates(1:end-13);
windowSize=5;
numWindows=length(breakDates)-windowSize;
rollingLogRets=cell(numWindows,1);
inRange=zeros(length(d),numWindows);
numDays=zeros(numWindows,1);
covMatWindow=cell(numWindows,1);
for i=1:numWindows
    %splitting into rolling windows
    inRange(:,i)=(d>breakDates(i))&(d<breakDates(i+windowSize));
    rows=find(inRange(:,i));
    rollingLogRets{i}=logRets(rows,:);
    covMatWindow{i}=cov(rollingLogRets{i});
    
    %finding MVP
    numDays(i)=length(rollingLogRets{i});
    H = 2.0 * covMatWindow{i};
    
    [weightsWindowShort{i}, sigmaWindowShort{i}] = quadprog(H, f, A_short, b_short, Aeq,beq, [], [], [], []);
    
    [weightsWindowNoShort{i}, sigmaWindowNoShort{i}] = quadprog(H, f, A_noShort, b_noShort, Aeq,beq, [], [], [], []);
end

