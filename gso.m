clc;
clear;
close all;

%% Problem Definition

model=CreateModel();

model.n=3;  % number of Handle Points

CostFunction=@(x) MyCost(x,model);    % Cost Function

nVar=model.n;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin.x=model.xmin;           % Lower Bound of Variables
VarMax.x=model.xmax;           % Upper Bound of Variables
VarMin.y=model.ymin;           % Lower Bound of Variables
VarMax.y=model.ymax;           % Upper Bound of Variables


%% GSO Parameters

MaxIt=10;          % Maximum Number of Iterations

nPop=150;           % Population Size (Swarm Size)

w=1;                % Inertia Weight
wdamp=0.98;         % Inertia Weight Damping Ratio
c1=1.5;             % Personal Learning Coefficient
c2=1.5;             % Global Learning Coefficient

%RANGE
range_init = 5.0;
range_boundary = 25.2;

%LUCIFERIN
luciferin_init = 5;
luciferin_decay = 0.4;
luciferin_enhancement = 0.6;

%Neighbors
k_neigh = 5;
beta = 0.005;

step_size = 5;

%% Initialization

% Create Empty Glowworm Structure
empty_glowworm.Position=[];
empty_glowworm.range=[];
empty_glowworm.luciferin=[];
empty_glowworm.Cost=[];
empty_glowworm.Sol=[];
empty_glowworm.neighbors=[];
empty_glowworm.Best.Position=[];
empty_glowworm.Best.Cost=[];
empty_glowworm.Best.Sol=[];

% Initialize Global Best
GlobalBest.Cost=inf;

% Create glowworms Matrix
glowworm=repmat(empty_glowworm,nPop,1);

% Initialization Loop
for i=1:nPop
	% Initialize Position
    if i > 1
        glowworm(i).Position=CreateRandomSolution(model);
    else
        % Straight line from source to destination
        xx = linspace(model.xs, model.xt, model.n+2);
        yy = linspace(model.ys, model.yt, model.n+2);
        glowworm(i).Position.x = xx(2:end-1);
        glowworm(i).Position.y = yy(2:end-1);
    end

    % Initialize luciferin
    glowworm(i).luciferin.x=repmat( luciferin_init , 1 , nVar);
    glowworm(i).luciferin.y=repmat( luciferin_init , 1 , nVar);
    
    %Initialize range 
    glowworm(i).range.x = repmat( range_init , 1 , nVar);
    glowworm(i).range.y = repmat( range_init , 1 , nVar);

    neighbors = [];
    % Evaluation
    [glowworm(i).Cost, glowworm(i).Sol]=CostFunction(glowworm(i).Position);
    
    % Update Personal Best
    glowworm(i).Best.Position=glowworm(i).Position;
    glowworm(i).Best.Cost=glowworm(i).Cost;
    glowworm(i).Best.Sol=glowworm(i).Sol;
    
    % Update Global Best
    if glowworm(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=glowworm(i).Best;
        
    end
end

% Array to Hold Best Cost Values at Each Iteration
BestCost=zeros(MaxIt,1);


for it=1:MaxIt
    
    for i=1:nPop
        
        % x Part

    	 % Update luciferin
        glowworm(i).luciferin.x = (1-luciferin_decay).*glowworm(i).luciferin.x + luciferin_enhancement.*glowworm(i).Cost;


        neighbors.x = [];
        for k =1:nPop

			dist = abs(glowworm(i).Position.x - glowworm(k).Position.x);
			%if it is in it's range of sigth and it's brightness is higher
            if all(dist ~= 0 ) & dist <= glowworm(i).range.x & glowworm(i).luciferin.x <= glowworm(k).luciferin.x
                neighbors.x = [neighbors.x ; k];
            end
        end
        
        if size(neighbors.x,2) > 0
         % find the node in the direction of which the glowworm should
         % follow
            li = glowworm(i).luciferin.x;
            sum_lk = sum(glowworm(i).luciferin.x);
            neighbors_index = size(neighbors.x,2);
            
            %calc probabilties for each neighbor been followed
            probs = zeros(1,neighbors_index);
            for j = 1:neighbors_index
                probs(j) = sum(glowworm(j).luciferin.x) - sum(li);
            end
            probs = probs./( sum_lk - (size(probs,2)*li));
            
            %calc prob range
            acc = 0;
            wheel = [];
            for val = 1:size(probs,2)
                acc = acc + probs(val);
                wheel = [wheel ; acc];
            end

            %wheel(-1) = 1 ;

            %randomly choice a value for wheel selection method
            rand_val = rand;
            following = 1;
            
            
            
            for k = 1:size(wheel,2)
                if rand_val <= wheel(k)
                    following = k;
                end
            end

            toward_index = following;
            
            %Position update 
            glowworms = glowworm(i).Position.x;
            toward = glowworm(toward_index).Position.x;
            
            normV = norm(toward - glowworms);
            if normV == 0 && isnan(normV)
                normV = 5; %step size 
            end
            
            new_position = glowworms + step_size.*(toward-glowworms)./normV;
            glowworm(i).Position.x = new_position;
        end
        
        for p = 1:nVar
             glowworm(i).range.x(p) = min(range_boundary,max(0,glowworm(i).range.x(p) + (beta*(k_neigh-size(neighbors.x,2)))));
        end
        
       
        
        
        
         % y Part

    	 % Update luciferin
        glowworm(i).luciferin.y = (1-luciferin_decay)*glowworm(i).luciferin.y + luciferin_enhancement*glowworm(i).Cost;

        %find neighbour 
        
        neighbors.y = [];
        for k =1:nPop

			dist = abs(glowworm(i).Position.y - glowworm(k).Position.y);
			%if it is in it's range of sigth and it's brightness is higher
            if all(dist ~= 0 ) & dist <= glowworm(i).range.x & glowworm(i).luciferin.x < glowworm(k).luciferin.x
                neighbors.y = [neighbors.y ; k];
            end
        end
        
        if size(neighbors.y,2) > 0
            
            % find the node in the direction of which the glowworm should follow
            li = glowworm(i).luciferin.y;
            sum_lk = sum(glowworm(i).luciferin.y);
            neighbors_index = size(neighbors.y,2);
            
            %calc probabilties for each neighbor been followed
            probs = zeros(1,neighbors_index);
            for j = 1:neighbors_index
                probs(j) = sum(glowworm(j).luciferin.y) - sum(li);
            end
            probs = probs./(sum_lk - (size(probs,2)*li));
            
            %calc prob range
            acc = 0;
            wheel = [];
            for val = 1:size(probs,2)
                acc = acc + probs(val);
                wheel = [wheel ; acc];
            end

            %wheel(-1) = 1 ;

            %randomly choice a value for wheel selection method
            rand_val = rand;
            following = 1;
            
            
            
            for k = 1:size(wheel,2)
                if rand_val <= wheel(k)
                    following = k;
                end
            end

            toward_index = following;
            
            %Position update 
            glowworms = glowworm(i).Position.y;
            toward = glowworm(toward_index).Position.y;
            
            normV = norm(toward-glowworms);
            if normV == 0 && isnan(normV)
                normV = 5; %step size 
            end
            
            new_position = glowworms + step_size.*(toward-glowworms)./normV;
            A = new_position;
            glowworm(i).Position.y = new_position;
            
        end
        for p = 1:nVar
             glowworm(i).range.y(p) = min(range_boundary,max(0,glowworm(i).range.y(p) + (beta*(k_neigh-size(neighbors.y,2)))));
        end
 
        
           % Evaluation
        [glowworm(i).Cost, glowworm(i).Sol]=CostFunction(glowworm(i).Position);
        
           % Update Personal Best
        if glowworm(i).Cost<glowworm(i).Best.Cost

            glowworm(i).Best.Position=glowworm(i).Position;
            glowworm(i).Best.Cost=glowworm(i).Cost;
            glowworm(i).Best.Sol=glowworm(i).Sol;
        end 
        % Update Global Best
        if glowworm(i).Best.Cost<GlobalBest.Cost
            GlobalBest=glowworm(i).Best;
        end
            
        
        
    end
    % Update Best Cost Ever Found
    BestCost(it)=GlobalBest.Cost;
    
     % Show Iteration Information
    if GlobalBest.Sol.IsFeasible
        Flag=' *';
    else
        Flag=[', Violation = ' num2str(GlobalBest.Sol.Violation)];
    end
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it)) Flag]);
    
    % Plot Solution
    figure(1);
    PlotSolution(GlobalBest.Sol,model);
    pause(0.01);
end

%% Results

figure;
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

        