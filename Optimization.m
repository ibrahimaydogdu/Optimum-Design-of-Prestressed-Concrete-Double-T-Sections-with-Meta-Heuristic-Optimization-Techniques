function [GlobalBest,GlobalParam,hisx,hisy] = Optimization(ub,lb,irun)
global Opt_Method
if Opt_Method == 1
    [GlobalBest,GlobalParam,hisx,hisy] = runABC(ub,lb,irun);
elseif Opt_Method == 2
     [GlobalBest, GlobalParam, hisx, hisy]= BSO(ub,lb,irun);
else
    [GlobalBest, GlobalParam, hisx, hisy]= SSO(ub,lb,irun);
end
end

%% ABC
function [GlobalBest,GlobalParam,hisx,hisy]= runABC(ub,lb,irun)
%global nu ffc g_conc
n=size(ub,1);

%% /* Control Parameters of ABC algorithm*/
NP=100; %/* The number of colony size (employed bees+onlooker bees)*/
FoodNumber=NP/2; %/*The number of food sources equals the half of the colony size*/
limit=50; %/*A food source which could not be improved through "limit" trials is abandoned by its employed bee*/
maxCycle=10000; %/*The number of cycles for foraging {a stopping criteria}*/
if irun==1
    disp('---------------------------------------------------------');
    disp('Artificial Bee Colony(ABC) Algorithm Parameters')
    fprintf('The number of food sources  = %i\n',FoodNumber);
    fprintf('Nectar limit of food source = %i\n',limit);
    fprintf('Maximum iteration Number    = %i\n',maxCycle);
    disp('---------------------------------------------------------');
end

%% First Part: Initialization

%reset trial counters
trial=zeros(1,FoodNumber);

% /*All food sources are initialized */
%/*Variables are initialized in the range [lb,ub]. If each parameter has different range, use arrays lb[j], ub[j] instead of lb and ub */
Foods=zeros(FoodNumber,n);
for i=1:FoodNumber
    for j=1:n
        Foods(i,j)=round(lb(j)+rand*(ub(j)-lb(j)));
    end
end

% Calculate nectar amount of the food source
Fitness=zeros(FoodNumber,1);ObjVal=zeros(FoodNumber,1);Const=zeros(FoodNumber,1);
FitBest=0;iter=0;hisx=zeros(1,1);hisy=zeros(1,1);icount=0;GlobalBest=1e40;
GlobalParam=zeros(n,1);
for i=1:FoodNumber
    x=Foods(i,:);
    [ObjVal(i),Const(i),Fitness(i),Foods(i,:)]= Computation(x);
    iter=iter+1;
    [icount,FitBest,GlobalBest,GlobalParam,hisx,hisy]= ...
        GlobalCheck(icount,Foods(i,:),Fitness(i),ObjVal(i),GlobalBest,...
        FitBest,Const(i),GlobalParam,iter,hisx,hisy);     
end

iter=1;
while ((iter <= maxCycle))

%% Second Part: EMPLOYED BEE PHASE
    for i=1:(FoodNumber)
        
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*n)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
        while(neighbour==i)
        	neighbour=fix(rand*(FoodNumber))+1;
        end
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=round(Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2);
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       if sol(Param2Change) > ub(Param2Change);sol(Param2Change) = ub(Param2Change);end
       if sol(Param2Change) < lb(Param2Change);sol(Param2Change) = lb(Param2Change);end 
        
        %evaluate new solution
        [ObjValSol,c,FitnessSol,sol]=Computation(sol);
        iter=iter+1;

       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
            Const(i)=c;
            
      %Global Check
           [icount,FitBest,GlobalBest,GlobalParam,hisx,hisy]= ...
           GlobalCheck(icount,sol,FitnessSol,ObjValSol,GlobalBest,FitBest,c,...
           GlobalParam,iter,hisx,hisy);
       else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end
               
    end

%% Mid Part: CalculateProbabilities
%/* A food source is chosen with the probability which is proportioal to its quality*/
%/*Different schemes can be used to calculate the probability values*/
%/*For example prob(i)=fitness(i)/sum(fitness)*/
%/*or in a way used in the metot below prob(i)=a*fitness(i)/max(fitness)+b*/
%/*probability values are calculated by using fitness values and normalized
%by dividing maximum fitness value*/
prob=(0.9.*Fitness./max(Fitness))+0.1;
  
%% Third Part: ONLOOKER BEE PHASE

i=1;
t=0;
while(t<FoodNumber)
    if(rand<prob(i))
        t=t+1;
        %/*The parameter to be changed is determined randomly*/
        Param2Change=fix(rand*n)+1;
        
        %/*A randomly chosen solution is used in producing a mutant solution of the solution i*/
        neighbour=fix(rand*(FoodNumber))+1;
       
        %/*Randomly selected solution must be different from the solution i*/        
        while(neighbour==i)
        	neighbour=fix(rand*(FoodNumber))+1;
        end
        
       sol=Foods(i,:);
       %  /*v_{ij}=x_{ij}+\phi_{ij}*(x_{kj}-x_{ij}) */
       sol(Param2Change)=round(Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2);
        
       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
       if sol(Param2Change) > ub(Param2Change);sol(Param2Change) = ub(Param2Change);end
       if sol(Param2Change) < lb(Param2Change);sol(Param2Change) = lb(Param2Change);end 
        
        %evaluate new solution
        [ObjValSol,c,FitnessSol,sol]=Computation(sol);
        iter=iter+1;
        
       % /*a greedy selection is applied between the current solution i and its mutant*/
       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
            Foods(i,:)=sol;
            Fitness(i)=FitnessSol;
            ObjVal(i)=ObjValSol;
            trial(i)=0;
            Const(i)=c;
            
      %Global Check
           [icount,FitBest,GlobalBest,GlobalParam,hisx,hisy]= ...
           GlobalCheck(icount,sol,FitnessSol,ObjValSol,GlobalBest,FitBest,c,...
           GlobalParam,iter,hisx,hisy);
       else
            trial(i)=trial(i)+1; %/*if the solution i can not be improved, increase its trial counter*/
       end
    end
    
    i=i+1;
    if (i==(FoodNumber)+1) 
        i=1;
    end 
end 
   
%% Fourth Part:SCOUT BEE PHASE 

%/*determine the food sources whose trial counter exceeds the "limit" value. 
%In Basic ABC, only one scout is allowed to occur in each cycle*/

ind=find(trial==max(trial));
ind=ind(end);
if (trial(ind)>limit)
    trial(ind)=0;
    for j=1:n
        sol(j)=round(lb(j)+rand*(ub(j)-lb(j)));
    end
    [ObjValSol,c,FitnessSol,sol]=Computation(sol);
    Foods(ind,:)=sol;
    Fitness(ind)=FitnessSol;
    ObjVal(ind)=ObjValSol;
    Const(ind)=c;
end

end % End of ABC
hisx(icount+1) = maxCycle;
hisy(icount+1) = GlobalBest;
if GlobalBest==1e40
    [~,i]=min(ObjVal);
    GlobalParam=Foods(i,:)';
end
end

%% BSO
function [GlobalBest, GlobalParam, hisx, hisy]= BSO(ub,lb,irun)
%warning('off','All warnings disableyed')
% n_p; population size
% n_d; number of dimension
% n_c: number of clusters
n_p=50;
n_d=size(ub,1);
n_c=5;
prob_one_cluster = 0.8; % probability for select one cluster to form new individual;
max_iteration=5000;
if irun==1
    disp('---------------------------------------------------------');
    disp('Brain Storming Algorithm Parameters')
    fprintf('Population Size=                %i\n',n_p);
    fprintf('Numberof Cluster=               %i\n',n_c);
    fprintf('probability for select one \n cluster to form new individual= %4.2f\n',prob_one_cluster);
    fprintf('Maximum iteration=              %i\n',max_iteration);
    disp('---------------------------------------------------------');
end
popu=zeros(n_p,n_d);popu_sorted=zeros(n_p,n_d);
Objective=zeros(n_p,1);Const=zeros(n_p,1);
FitBest=0;hisx=zeros(1,1);hisy=zeros(1,1);icount=0;GlobalBest=1e40;
GlobalParam=zeros(n_d,1);
fitness_popu = 1000000*ones(n_p,1);  % store fitness value for each individual
fitness_popu_sorted = 1000000*ones(n_p,1);  % store  fitness value for each sorted individual
n_iteration = 0; % current iteration number
for i=1:n_p
    for j=1:n_d
        popu(i,j)=round(lb(j) + rand*(ub(j) - lb(j))); % initialize the population of individuals
        popu_sorted(i,j)  =round(lb(j) + rand*(ub(j) - lb(j))); % initialize the  population of individuals sorted according to clusters  
    end
    
    % calculate fitness for each individual in the initialized population
    [Objective(i), Const(i), fitness_popu(i),popu(i,:)]=Computation(popu(i,:));
    n_iteration=n_iteration+1;
    [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]= ...
    GlobalCheck(icount,popu(i,:),fitness_popu(i),Objective(i),GlobalBest,...
    FitBest,Const(i),GlobalParam,n_iteration,hisx,hisy); 
end


% initialize cluster probability to be zeros
prob = zeros(n_c,1);
best = zeros(n_c,1);  % index of best individual in each cluster

centers=zeros(n_c,n_d);centers_copy=zeros(n_c,n_d);
for i=1:n_c
    for j=1:n_d
        centers(i,j)=round(lb(j) + rand*(ub(j) - lb(j))); % initialize best individual in each cluster
        centers_copy(i,j)  =round(lb(j) + rand*(ub(j) - lb(j)));% initialize best individual-COPY in each cluster FOR the purpose of introduce random best
    end
end

indi_temp = zeros(1,n_d);  % store temperary individual




while n_iteration < max_iteration
    
    cluster=kmeans(popu, n_c,'Distance','cityblock','start',centers,'emptyaction','singleton'); % k-mean cluster
    
    % clustering    
    fit_values = 100000000000000000000000000.0*ones(n_c,1);  % assign a initial big fitness value  as best fitness for each cluster in minimization problems

    number_in_cluster = zeros(n_c,1);  % initialize 0 individual in each cluster
           
    for idx = 1:n_p
        number_in_cluster(cluster(idx,1),1)= number_in_cluster(cluster(idx,1),1) + 1;
               
        % find the best individual in each cluster
        if fit_values(cluster(idx,1),1) > fitness_popu(idx,1)  % minimization
            fit_values(cluster(idx,1),1) = fitness_popu(idx,1);
            best(cluster(idx,1),1) = idx;
        end
            
    end
    
    % form population sorted according to clusters
    counter_cluster = zeros(n_c,1);  % initialize cluster counter to be 0 
    
    acculate_num_cluster = zeros(n_c,1);  % initialize accumulated number of individuals in previous clusters
    
    for idx =2:n_c
        acculate_num_cluster(idx,1) = acculate_num_cluster((idx-1),1) + number_in_cluster((idx-1),1);
    end
    
    
    %start form sorted population
    for idx = 1:n_p
        counter_cluster(cluster(idx,1),1) = counter_cluster(cluster(idx,1),1) + 1 ;
        temIdx = acculate_num_cluster(cluster(idx,1),1) +  counter_cluster(cluster(idx,1),1);
        popu_sorted(temIdx,:) = popu(idx,:);
        fitness_popu_sorted(temIdx,1) = fitness_popu(idx,1);
    end
       
    % record the best individual in each cluster
    for idx = 1:n_c
        centers(idx,:) = popu(best(idx,1),:);        
    end
    
    centers_copy = centers;  % make a copy
    
    if (rand() < 0.2) %  select one cluster center to be replaced by a randomly generated center
        cenIdx = ceil(rand()*n_c);
        for j=1:n_d
            centers(cenIdx,j) = round(lb(j) + (ub(j) - lb(j)) * rand);
        end
    end 
           
    % calculate cluster probabilities based on number of individuals in
    % each cluster
    for idx = 1:n_c
        prob(idx,1) = number_in_cluster(idx,1)/n_p;
        if idx > 1
            prob(idx,1) = prob(idx,1) + prob(idx-1,1);
        end
    end
    
    % generate n_p new individuals by adding Gaussian random values
                   
    for idx = 1:n_p
        for j=1:n_d
        r_1 = rand();  % probability for select one cluster to form new individual
        if r_1 < prob_one_cluster % select one cluster
            r = rand();
            for idj = 1:n_c
                if r < prob(idj,1)                      
                    if rand() < 0.4  % use the center
                       indi_temp(1,j) = centers(idj,j); 
                    else % use one randomly selected  cluster
                        indi_1 = acculate_num_cluster(idj,1) + ceil(rand() * number_in_cluster(idj,1));
                        indi_temp(1,j) = popu_sorted(indi_1,j);  
                    end
                    break
                end
            end
        else % select two clusters
            % pick two clusters 
            cluster_1 = ceil(rand() * n_c);
            indi_1 = acculate_num_cluster(cluster_1,1) + ceil(rand() * number_in_cluster(cluster_1,1));
            
            cluster_2 = ceil(rand() * n_c);
            indi_2 = acculate_num_cluster(cluster_2,1) + ceil(rand() * number_in_cluster(cluster_2,1));
            
            tem = rand();
            if rand() < 0.5 %use center
                indi_temp(1,j) = tem * centers(cluster_1,j) + (1-tem) * centers(cluster_2,j); 
            else   % use randomly selected individuals from each cluster            
                indi_temp(1,j) = tem * popu_sorted(indi_1,j) + (1-tem) * popu_sorted(indi_2,j); 
            end
            
            
        end
        end
        %stepSize = ones(1,n_d); % effecting the step size of generating
        %new individuals by adding random values kaldýrýldý 26.01.17
%         stepSize = logsig((0.5*max_iteration - n_iteration)/(20*n_p)) * rand; %26.01.17 this part is updated standart stepsize causes divergence 
%         indi_temp(1,dimen) = indi_temp(1,dimen) + stepSize* normrnd(0,1,1)*(rang_r(dimen)-rang_l(dimen));
        if rand>fitness_popu(idx,1)/max(fitness_popu)
             dimen=fix(rand*n_d)+1;
             indi_temp(1,dimen)=round(lb(dimen) + rand*(ub(dimen) - lb(dimen)));
        end

        %Upper lower boundry check
        for j=1:n_d
            if indi_temp(1,j) > ub(j);indi_temp(1,j) = ub(j);end
            if indi_temp(1,j) < lb(j);indi_temp(1,j) = lb(j);end
            indi_temp(1,j)=round(indi_temp(1,j));
        end
        % if better than the previous one, replace it
        [Obj, c, fv,indi_temp(1,:)]=Computation(indi_temp(1,:));
        n_iteration = n_iteration +1;
        if fv > fitness_popu(idx,1)  % better than the previous one, replace
            fitness_popu(idx,1) = fv;
            popu(idx,:) = indi_temp(1,:);
            Objective(idx)=Obj;
            Const(idx)=c;        
            [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]= ...
            GlobalCheck(icount,popu(idx,:),fitness_popu(idx,1),Objective(idx),GlobalBest,...
            FitBest,Const(idx),GlobalParam,n_iteration,hisx,hisy); 
        end
        
    end
    % Clean same solutions in the colony
    [popu, Const, fitness_popu, Objective, GlobalBest, FitBest, GlobalParam,...
    n_iteration, hisx, hisy]= Clear_Dump(popu,Const,fitness_popu,Objective,...
    GlobalBest,FitBest,GlobalParam,n_iteration,hisx,hisy,icount,ub,lb);  

    % keep the best for each cluster
    for idx = 1:n_c
        popu(best(idx,1),:) = centers_copy(idx,:);  
        fitness_popu(best(idx,1),1) = fit_values(idx,1);
    end
end
hisx(icount+1) = max_iteration;
hisy(icount+1) = GlobalBest;
end


%% SSO
function [GlobalBest, GlobalParam, hisx, hisy]= SSO(ub,lb,irun)
%Social Spider Algorithm Matlab Code
%The code written by Dr. Ibrahim Aydogdu 18.06.2015
n=size(ub,1);
%Algorithm Parameters
%Number of spider :NofSpider
NofSpider=50;
%Probability of selection: Pf
Pf=0.3;
%Maximum iteration number: MaxIter
MaxIter=5000;
if irun==1
    disp('---------------------------------------------------------');
    disp('Social Spider Algorithm Parameters')
    fprintf('Number of spider=         %i\n',NofSpider);
    fprintf('Probability of selection= %4.2f\n',Pf);
    fprintf('Maximum iteration=        %4.2f\n',MaxIter);
    disp('---------------------------------------------------------');
end
%% First Part: Initialization
%Define number of male female spiders
 Nf = int8((0.9 - rand * 0.25) * NofSpider);
 Nm = NofSpider - Nf;
 
 %Determine initial positions of spiders randomly
Spider=zeros(NofSpider,n);
for i=1:NofSpider
    for j=1:n
        Spider(i,j)=round(lb(j)+rand*(ub(j)-lb(j)));
    end
end

%Calculating mating radius and maximum distance
 r=zeros(n);
 for i=1:n
    r(i)=(ub(i)-lb(i))/2;
 end
 distmax=sqrt(sum((ub-lb).^2));

%Measurement of performance of spider
Fitness=zeros(NofSpider,1);Objective=zeros(NofSpider,1);Const=zeros(NofSpider,1);
FitBest=0;iter=0;hisx=zeros(1,1);hisy=zeros(1,1);icount=0;GlobalBest=1e40;
GlobalParam=zeros(n,1);
for i=1:NofSpider
        [Objective(i), Const(i), Fitness(i),Spider(i,:)]=Computation(Spider(i,:));
        iter=iter+1;
        [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]= ...
        GlobalCheck(icount,Spider(i,:),Fitness(i),Objective(i),GlobalBest,...
        FitBest,Const(i),GlobalParam,iter,hisx,hisy);     
end


 while iter<MaxIter
% Generate f and m matrices
      f=zeros(Nf,n);m=zeros(Nm,n);
      for i=1:NofSpider
          if i<=Nf
              f(i,:)=Spider(i,:);
          else
              m(i-Nf,:)=Spider(i,:);
          end
      end
      
%Saving worst and best spider
    [BestS, iBestS]=max(Fitness);
    [WorstS, iWorstS]=min(Fitness);

% Calculate weight of Spiders
    w=zeros(NofSpider,1);
    for i=1:NofSpider
        w(i)=(WorstS-Fitness(i))/(WorstS-BestS);
    end
%% Second Part: Calculation of vibration
%female
    vibc=zeros(Nf,1);bc=zeros(Nf,1);vibb=zeros(Nf,1);Dist=zeros(Nf,1);
    for i = 1:Nf
    %Closest Female                
        Dist(:)=1e39;
        for k = 1:Nf
        	if i ~= k;Dist(k)=norm(f(i,:)-f(k,:))/distmax;end
        end
        [Distmin, bc(i)]=min(Dist);
        if Fitness(i) > Fitness(bc(i))
        	vibc(i) = w(i) * exp(-Distmin^2);
        else
            vibc(i) = w(bc(i)) * exp(-Distmin^2);
        end     
        if Fitness(i)~= BestS
    %Best female
            Dist=norm(f(i,:)-Spider(iBestS,:))/distmax; 
            vibb(i) = w(i) * exp(-Dist * Dist);
        end
    end
    
%Male spiders
    vibf=zeros(Nm,1);bf=zeros(Nm,1);
    for i = 1:Nm
            Dist(:)=1e39;   
            for k = 1:Nf
                Dist(k)=norm(m(i,:)-f(k,:))/distmax;
            end
            [Distmin, bf(i)]=min(Dist);
            vibf(i) = w(bf(i)) * exp(-Distmin * Distmin);
    end
 %% Third Part:Movement Operators
 %Female cooperative operator
	for i = 1:Nf
        if Fitness(i)~= BestS
            for j = 1:n
                if rand < Pf 
                	f(i, j) = f(i, j) + round(rand * vibc(i) * (f(bc(i), j) - f(i, j)) + rand * vibb(i)*(Spider(iBestS, j) - f(i, j)) + rand * (rand - 0.5));
                else
                    f(i, j) = f(i, j) - round(rand * vibc(i) * (f(bc(i), j) - f(i, j)) - rand * vibb(i)*(Spider(iBestS, j) - f(i, j)) + rand * (rand - 0.5)); 

                end
                if f(i, j) > ub(j);f(i, j) = ub(j);end
                if f(i, j) < lb(j);f(i, j) = lb(j);end        
            end
        end
	end
    
%  Male cooperative operator
    % Sorting
        WM=zeros(Nm,1);
        for i=1:Nm;WM(i)=w(Nf+i);end
        SWM=sort(WM,'descend');
        WNfm = SWM(round((Nm + 1) / 2));
        SmWNfh=zeros(n,1);
        for j=1:n
        	SmWNfh(j)=sum(m(:, j).* WM);
        end
        SWNfh=sum(WM);
        for i = 1:Nm
             for j = 1:n
                 if WM(i) > WNfm
                     m(i, j) = m(i, j) + round(rand * vibf(i) * (f(bf(i), j) - m(i, j)) + rand * (rand - 0.5));
                 else
                     m(i, j) = m(i, j) + round(rand * (SmWNfh(j) / SWNfh -m(i, j)));                    
                 end
                 if m(i, j) > ub(j);m(i, j) = ub(j);end
                 if m(i, j) < lb(j);m(i, j) = lb(j);end
             end
        end
        iter=iter+1;
        etkim=zeros(Nm, n, Nf);Ps=zeros(Nf + 1); %SPs=zeros(Nf + 1);ind=zeros(Nf+1);
        NewSp=zeros(n,1);
        
% After movement update spider matrix
        Spider=[f;m];
%% Fourth Part:Mating operator
         for i = 1:Nm
             if WM(i) > WNfm
                 for j = 1:n
                    ik = 0;
                    for k = 1:Nf
                        if f(k, j) <= m(i, j) + r(j) && f(k, j) >= m(i, j) - r(j)
                            ik = ik + 1;
                            etkim(i, j, ik) = k;
                        end
                    end
                     if ik > 0
% Calculation of probability
                         TPs = 0;
                        for ik1 = 1:ik
                            TPs = TPs + w(etkim(i, j, ik1));
                        end
                        for ik1 = 1:ik
                            Ps(ik1) = w(etkim(i, j, ik1)) / TPs;
                        end
                        Ps(ik + 1) = WM(i) / TPs;
% Roulette method
                        [SPs, ind]=sort(Ps,'descend');
                        ik1=1;
                         TPs = SPs(ik1);
                         while  TPs<rand
                             ik1=ik1+1;
                             TPs = TPs + SPs(ik1);
                         end
% Generate New_spider
                         if ind(ik1) > ik
                         	NewSp(j) = m(i, j);
                         else
                         	NewSp(j) = f(etkim(i, j, ind(ik1)), j);
                         end
                     else
  
                        NewSp(j) = Spider(iWorstS, j);
                     end
                    if NewSp(j) > ub(j);NewSp(j) = ub(j);end
                    if NewSp(j) < lb(j);NewSp(j) = lb(j);end
                 end
                 
% Check performance of new spider and replacement procedure
                [Obj, c, Fit,NewSp]=Computation(NewSp);
                iter=iter+1;
                 if Fit > WorstS
                    Spider(iWorstS, :) = NewSp(:);
                    Fitness(iWorstS)=Fit;
                    Objective(iWorstS)=Obj;
                    Const(iWorstS)=c;
                 [icount, FitBest, GlobalBest, GlobalParam,  hisx, hisy]= ...
                 GlobalCheck(icount,NewSp,Fit,Obj,GlobalBest,FitBest,c,...
                 GlobalParam,iter,hisx,hisy);
                 end
             end
         end
%% Fifth Part: Update Spider colony after movement and mating
% Clean same solutions in the colony
    [Spider, Const, Fitness, Objective, GlobalBest, FitBest, GlobalParam,...
    iter, hisx, hisy]= Clear_Dump(Spider,Const,Fitness,Objective,...
    GlobalBest,FitBest,GlobalParam,iter,hisx,hisy,icount,ub,lb);
% Update Colony
    for i=1:NofSpider
        [Objective(i), Const(i), Fitness(i),Spider(i,:)]=Computation(Spider(i,:));
        iter=iter+1;
        [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]= ...
        GlobalCheck(icount,Spider(i,:),Fitness(i),Objective(i),GlobalBest,...
        FitBest,Const(i),GlobalParam,iter,hisx,hisy);     
    end       
 end
 hisx(icount+1) = MaxIter;
 hisy(icount+1) = GlobalBest;
end


%% Global Check
function [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]= GlobalCheck(icount,x,Fit,Obj,GlobalBest,FitBest,c,GlobalParam,iter,hisx,hisy)
if Fit>FitBest  && c<=0
    icount=icount+1;
    GlobalBest=Obj;
    GlobalParam(:)=x(:);
    FitBest=Fit;
    fprintf('iter:%i GlobalBest:%6.3f TotalPenaly:%6.3f\n', iter,  GlobalBest, c)
    hisx(icount)=iter;hisy(icount)=GlobalBest;
end
end


%% Clear Dump
function [Spider, Const, Fitness, Objective, GlobalBest, FitBest, GlobalParam,iter, hisx, hisy] = Clear_Dump(Spider,Const,Fitness,Objective,GlobalBest,FitBest,GlobalParam,iter,hisx,hisy,icount,ub,lb)
n=size(Spider,2);NofSpider=size(Objective,1);
for i=1:NofSpider
    for k=1:NofSpider
        if i~=k && Objective(i)==Objective(k) %iki ayný tasarým þartý
            for j=1:n
                Spider(k,j)=round(lb(j)+rand*(ub(j)-lb(j)));     
            end
            [Objective(k), Const(k), Fitness(k),Spider(k,:)]=Computation(Spider(k,:));
            iter=iter+1; 
            [icount, FitBest, GlobalBest, GlobalParam, hisx, hisy]=...
                GlobalCheck(icount,Spider(k,:),Fitness(k),Objective(k),GlobalBest,...
                FitBest,Const(k),GlobalParam,iter,hisx,hisy);
        end
    end
end

end
