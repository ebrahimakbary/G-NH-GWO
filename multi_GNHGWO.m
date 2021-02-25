clc
clear

global NFE
NFE=0;

nPop=30;    % Number of search agents (Population Number)
MaxIt=1000; % Maximum number of iterations
nVar=30;    % Number of Optimization Variables

Total_Runs=25;

% Pre-allocating vectors and matrices
Cost_Rsult=nan(1,Total_Runs);
Rsult=nan(Total_Runs,MaxIt);
Mean=nan(1,14);
Best=Mean;
Std=Mean;
nfe=Mean;
fitness1=nan(1,nPop);

for nFun=1:14
    NFE=0;
    CostFunction=@(x,nFun) Cost(x,nFun);        % Cost Function
    for run_no=1:Total_Runs
        
        %% Problem Definition
        
        VarMin=-100;             % Decision Variables Lower Bound
        if nFun==7
            VarMin=-600;             % Decision Variables Lower Bound
        end
        if nFun==8
            VarMin=-32;             % Decision Variables Lower Bound
        end
        if nFun==9
            VarMin=-5;             % Decision Variables Lower Bound
        end
        if nFun==10
            VarMin=-5;             % Decision Variables Lower Bound
        end
        if nFun==11
            VarMin=-0.5;             % Decision Variables Lower Bound
        end
        if nFun==12
            VarMin=-pi;             % Decision Variables Lower Bound
        end
        if nFun==14
            VarMin=-100;             % Decision Variables Lower Bound
        end
        VarMax= -VarMin;             % Decision Variables Upper Bound
        if nFun==13
            VarMin=-3;             % Decision Variables Lower Bound
            VarMax= 1;             % Decision Variables Upper Bound
        end
        
        %%   Grey Wold Optimizer (GWO)

        % Initialize Best Solution (Alpha) which will be used for archiving
        Alpha_pos=zeros(1,nVar);
        Alpha_score=inf;
                
        % Initialize the positions of search agents
        Positions=rand(nPop,nVar).*(VarMax-VarMin)+VarMin;
        Positions1=rand(nPop,nVar).*(VarMax-VarMin)+VarMin;
        BestCosts=zeros(1,MaxIt);

        fitness(1:nPop)=inf;

        iter=0;  % Loop counter
        
        %% Main loop
        while iter<MaxIt
            for i=1:nPop
                
                % Return back the search agents that go beyond the boundaries of the search space
                Flag4ub=Positions1(i,:)>VarMax;
                Flag4lb=Positions1(i,:)<VarMin;
                Positions1(i,:)=(Positions1(i,:).*(~(Flag4ub+Flag4lb)))+VarMax.*Flag4ub+VarMin.*Flag4lb;
                
                % Calculate objective function for each search agent
                fitness1(i)= CostFunction(Positions1(i,:), nFun);
                
                %  Grey Wolves
                if fitness1(i)<fitness(i)
                    Positions(i,:)=Positions1(i,:);
                    fitness(i) =fitness1(i) ;
                end
                
                % Update Best Solution (Alpha) for archiving
                if fitness(i)<Alpha_score
                    Alpha_score=fitness(i);
                    Alpha_pos=Positions(i,:);
                end
                
            end
            
            a=2-(iter*((2)/MaxIt));  % a decreases linearly fron 2 to 0
            
            % Update the Position of all search agents
            for i=1:nPop
                for j=1:nVar
                    
                    GGG=randperm(nPop-1,3);
                    ind1= GGG>=i;
                    GGG(ind1)=GGG(ind1)+1;
                    m1=GGG(1);
                    m2=GGG(2);
                    m3=GGG(3);
                    
                    r1=rand;
                    r2=rand;
                    
                    A1=2*a*r1-a;
                    C1=2*r2;
                    
                    D_alpha=abs(C1*Positions(m1,j)-Positions(i,j));
                    X1=Positions(m1,j)-A1*D_alpha;
                    
                    r1=rand;
                    r2=rand;
                    
                    A2=2*a*r1-a;
                    C2=2*r2;
                    
                    D_beta=abs(C2*Positions(m2,j)-Positions(i,j));
                    X2=Positions(m2,j)-A2*D_beta;
                    
                    r1=rand;
                    r2=rand;
                    
                    A3=2*a*r1-a;
                    C3=2*r2;
                    
                    D_delta=abs(C3*Positions(m3,j)-Positions(i,j));
                    X3=Positions(m3,j)-A3*D_delta;
                    Positions1(i,j)=(X1+X2+X3)/3;
                    
                end
            end
            
            iter=iter+1;
            BestCosts(iter)=Alpha_score;
            
            fprintf('Func No= %-2.0f,  Run No= %-2.0f,  Iter= %g, Best Cost = %g\n',nFun,run_no,iter,Alpha_score);
            
        end
        
        %%% Results
        
        Cost_Rsult(1, run_no)=Alpha_score;
        Rsult(run_no,:)= BestCosts;
    end
    
    Mean(nFun)=mean(Cost_Rsult);
    Best(nFun)=min(Cost_Rsult);
    Std(nFun)=std(Cost_Rsult);
    nfe(nFun)=NFE;
    
    disp('==================================================')
    fprintf('Runned up to fun %g \n\n',nFun)
    fprintf('%6s%8s\n','nFun','Mean')
    for i=1:nFun
        fprintf('%4.0f     %g\n',i,Mean(i))
    end
    fprintf('\n');  
    disp('==================================================')
    fprintf('\n');  

    hold on
    plot(log(mean(Rsult)),'k','LineWidth',4); hold on

end

save Total_Results_GNHGWO
