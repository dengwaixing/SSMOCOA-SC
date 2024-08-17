function [ps,pf]=SSMMCOASC(func_name,VRmin,VRmax,n_obj,Particle_Number,maxgen,Q,P,c1)

% Dimension: n_var --- dimensions of decision space
%            n_obj --- dimensions of objective space
%% Input:
%                      Dimension                    Description
%      func_name       1 x length function name     the name of test function     
%      VRmin           1 x n_var                    low bound of decision variable
%      VRmax           1 x n_var                    up bound of decision variable
%      n_obj           1 x 1                        dimensions of objective space
%      Particle_Number 1 x 1                        population size
%      Max_Gen         1 x 1                        maximum  generations

%% Output:
%                     Description
%      ps             Pareto set
%      pf             Pareto front


%% Initialize parameters

sigma=0.1;
n_PBA=5;
n_var=size(VRmin,2);
%Canonical point set initialisation. 
tmp1 = [1: Particle_Number]'*ones(1,n_var );
Ind = [1: n_var];
prime1 = primes(100*n_var);
[p,q]=find(prime1 >= (2*n_var+3));
tmp2 = (2*pi.*Ind)/prime1(1,q(1));
tmp2 = 2*cos(tmp2);
tmp2 = ones(Particle_Number,1)*tmp2;
GD = tmp1.*tmp2;
GD = mod(GD,1);
GD=VRmin+(VRmax-VRmin).*GD;
plot(GD(:,1),GD(:,2),'*');
pos=GD;
%% Evaluate the population
  
 fitness=zeros(Particle_Number,n_obj);
 for i=1:Particle_Number
      fitness(i,:)=feval(func_name,pos(i,:));
 end          
particle=[pos,fitness];    
row_of_cell=ones(1,Particle_Number); 
col_of_cell=size(particle,2);       
PBA=mat2cell(particle,row_of_cell,col_of_cell);%Archive A

for gcount=1:maxgen
    num_clusters=round( Particle_Number/Q);
    disp(gcount)    
    fitpop=non_domination_combine_cd_sort1(particle(:,1:n_var+n_obj), n_obj, n_var,c1);
    spop=[];
    W = exp(-squareform(pdist(fitpop(:,1:n_var))).^2/(2*sigma^2));
    C = spectral1(W,num_clusters);%Spectral clustering methods to form subpopulations
    for j=1:num_clusters  
        K=find(C==j);
        spop(j).pop=fitpop(K,:);
        POP=non_domination_combine_cd_sort1(spop(j).pop(:,1:n_var+n_obj), n_obj, n_var,c1);
        spop(j).species=POP(1,1:n_var);
    end
    rank1=fitpop( fitpop(:,n_obj + n_var + 1)==1,:);
    rank2=fitpop( fitpop(:,n_obj + n_var + 1)==2,:);
    if size(rank2,1)<=10
             xelist=[rank1(1:size(rank1,1),:);rank2];
    else
               xelist=[rank1(1:size(rank1,1),:);rank2(1:round(size(rank2,1)/3),:)];
    end
    for i =1:num_clusters 
        for t=1:size(spop(i).pop,1)
            for j=1:size(particle,1)
                if(spop(i).pop(t,1:n_var)==particle(j,1:n_var))
                    v=j;
                end
            end
       
         I=rand()*1+P;
        c=(maxgen-gcount)/maxgen;
        if t<=size(spop(i).pop,1)*1/2
             iguana= spop(i).species;
             if gcount<=maxgen*1/5
                   dis= fitpop(1:size(xelist,1),1:n_var+n_obj);
              else
                 dis= fitpop(1:size(rank1,1),1:n_var+n_obj);
              end
               if size(rank1,1)<5
                    dis=dis(1:size(rank1,1),:);
               else
                   dis=dis(1:5,:);
                   
               end
              poptemp=repmat( spop(i).pop(t,1:n_var),size(  dis,1),1);
              distance=sqrt(sum((dis(:,1:n_var)-poptemp(:,1:n_var)).*(dis(:,1:n_var)-poptemp(:,1:n_var)),2));% popsize*1
             [~,ind_dis]=sort(distance');% sorting the distance to p_i th particle in desend order
              dis1=fitpop(ind_dis(1,1),1:n_var);
             spop(i).pop(t,1:n_var)=spop(i).pop(t,1:n_var)+c*rand(1,n_var).*(iguana-spop(i).pop(t,1:n_var))+(1-c)*rand(1,n_var).*(I*dis1-spop(i).pop(t,1:n_var));   
        
        else
            if gcount<=maxgen*1/5
                
              
                  iguana=(VRmax+0.3)-rand(1,n_var).*(VRmax-VRmin);
                  iguana=max( iguana,VRmin); 
                  iguana= min(iguana,VRmax);
            else
          xelist=fitpop(1:size(rank1,1),1:n_var+n_obj);
          xelist=non_domination_combine_cd_sort1(xelist(:,1:n_var+n_obj), n_obj, n_var,c1);
          xelist= xelist(1:2,:);
          poptemp=repmat( spop(i).pop(t,1:n_var),size(xelist ,1),1);
          distance=sqrt(sum((xelist(:,1:n_var)-poptemp(:,1:n_var)).*(xelist(:,1:n_var)-poptemp(:,1:n_var)),2));% popsize*1
          [~,ind_dis]=sort(distance');% sorting the distance to p_i th particle in desend order
          iguana=fitpop(ind_dis(1,1),1:n_var);
      
           end
           L_fitness=reshape(feval(func_name,iguana),1,n_obj); 
    
         if n_obj<3
          
            if particle(v,n_var+1)<=L_fitness(1,1)&&particle(v,n_var+2)<=L_fitness(1,2)      
                  spop(i).pop(t,1:n_var)=spop(i).pop(t,1:n_var)+(1-c)*rand(1,n_var).*(I*spop(i).pop(t,1:n_var)-iguana); % Eq. (6)
           
            else
                  spop(i).pop(t,1:n_var)=spop(i).pop(t,1:n_var)+(1-c)*rand(1,n_var).*(I*iguana-spop(i).pop(t,1:n_var));    
            end
         else
                 if particle(v,n_var+1)<=L_fitness(1,1)&&particle(v,n_var+2)<=L_fitness(1,2)&&particle(v,n_var+3)<=L_fitness(1,3)      
                  spop(i).pop(t,1:n_var)=spop(i).pop(t,1:n_var)+rand()*(1-c)*(spop(i).pop(t,1:n_var)-iguana); % Eq. (6)
           
            else
                  spop(i).pop(t,1:n_var)=spop(i).pop(t,1:n_var)+rand()*(1-c)*(iguana-I.*spop(i).pop(t,1:n_var));    
            end
        end
             spop(i).pop(t,1:n_var) = max(  spop(i).pop(t,1:n_var),VRmin); 
             spop(i).pop(t,1:n_var) = min(  spop(i).pop(t,1:n_var),VRmax);
             spop(i).new_fitness(t,1:n_obj)=feval(func_name,spop(i).pop(t,1:n_var));
          
             if n_obj<3
          if spop(i).pop(t,n_var+n_obj+1)==1||spop(i).pop(t,n_var+n_obj+1)==2
          if spop(i).new_fitness(t,1)<=particle(v,n_var+1)&&spop(i).new_fitness(t,2)<=particle(v,n_var+2)
              particle(v,:)=[spop(i).pop(t,1:n_var),spop(i).new_fitness(t,1:n_obj)];
             
          end
          else
             particle(v,:)=[spop(i).pop(t,1:n_var),spop(i).new_fitness(t,1:n_obj)];
          end
             else
               if spop(i).pop(t,n_var+n_obj+1)==1||spop(i).pop(t,n_var+n_obj+1)==2
                if particle(v,n_var+1)>=L_fitness(1,1)&&particle(v,n_var+2)>=L_fitness(1,2)&&particle(v,n_var+3)>=L_fitness(1,3)     
                    particle(v,:)=[spop(i).pop(t,1:n_var),spop(i).new_fitness(t,1:n_obj)];
                 
                end
               else
                    particle(v,:)=[spop(i).pop(t,1:n_var),spop(i).new_fitness(t,1:n_obj)];
          end
                   
          end
              
  
             PBA_v=PBA{v,1};  
           PBA_v=[PBA_v(:,1:n_var+n_obj);[spop(i).pop(t,1:n_var),spop(i).new_fitness(t,1:n_obj)]];
            PBA_v=unique(PBA_v,'row');
           PBA_v= non_domination_combine_cd_sort1(PBA_v(:,1:n_var+n_obj), n_obj, n_var,c1);
           if size(PBA_v,1)>n_PBA
                 PBA{v,1}=PBA_v(1:n_PBA,:);
           else
                 PBA{v,1}=PBA_v;
           end
        end
        end
    end
     for i=1:Particle_Number
         LO_LOCAL=VRmin/gcount;% Eq(9)
        HI_LOCAL=VRmax/gcount;% Eq(10)
        new_particle1(i,1:n_var)=particle(i,1:n_var)+(1-2*rand(1,n_var)).* (LO_LOCAL+rand(1,n_var) .* (HI_LOCAL-LO_LOCAL)); % Eq. (8)
        new_particle1(i,1:n_var) = max( new_particle1(i,1:n_var),VRmin); new_particle1(i,1:n_var) = min( new_particle1(i,1:n_var),VRmax);
        new_fitness1(i,1:n_obj)=feval(func_name,new_particle1(i,1:n_var));
      if n_obj<3  
          if new_fitness1(i,1)<=particle(i,n_var+1)&new_fitness1(i,2)<=particle(i,n_var+2)
              particle(i,:)=[new_particle1(i,1:n_var),new_fitness1(i,1:n_obj)];
          end
          else
                if particle(i,n_var+1)>=new_fitness1(1,1)&&particle(i,n_var+2)>=new_fitness1(1,2)&&particle(i,n_var+3)>=new_fitness1(1,3)     
                    particle(i,:)=[new_particle1(i,1:n_var),new_fitness1(i,1:n_obj)];
                end
    end
    
           PBA_v=PBA{i,1};  
           PBA_v=[PBA_v(:,1:n_var+n_obj);[new_particle1(i,1:n_var),new_fitness1(i,1:n_obj)]];
             PBA_v=unique(PBA_v,'row');
           PBA_v= non_domination_combine_cd_sort1(PBA_v(:,1:n_var+n_obj), n_obj, n_var,c1);
           if size(PBA_v,1)>n_PBA
                 PBA{i,1}=PBA_v(1:n_PBA,:);
           else
                 PBA{i,1}=PBA_v;
           end
     end
end
  tempEXA=cell2mat(PBA); 
    tempEXA=non_domination_combine_cd_sort1(tempEXA(:,1:n_var+n_obj),n_obj,n_var,c1);
     if size(tempEXA,1)>Particle_Number
         EXA=tempEXA(1:Particle_Number,:);
     else
         EXA=tempEXA;
     end
  
   ps=EXA(:,1:n_var);
   pf=EXA(:,n_var+1:n_var+n_obj);

end
    




