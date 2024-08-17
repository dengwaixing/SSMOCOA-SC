function [population]=non_domination_combine_cd_sort(x,n_obj,n_var,c1)
% computing crowding distance of new possible next generation patents considering already selected parents for next generation
% Input  x            the points need to be ranked: popsize*(n_var+n_obj)
%        n_obj        number of objective 
%        n_var        number of varibles
%        Rate_F1      the rate selected particles in the first rank
%        Num_sele     the number of particles need to be selected from x
% Output sorted_based_on_front      the rank value is added to the
%     P_selected      the selected population

    [N_particle, ~] = size(x);% Obtain the number of particles
    front = 1;% Initialize the front number to 1.
    F(front).f = [];
    individual = [];
    population=[];
%% Non-Dominated sort. 

    for i = 1 : N_particle
        % Number of individuals that dominate this individual
        individual(i).n = 0; 
        % Individuals which this individual dominate
        individual(i).p = [];
        for j = 1 : N_particle
            dom_less = 0;
            dom_equal = 0;
            dom_more = 0;
            for k = 1 : n_obj
                if (x(i,n_var + k) < x(j,n_var + k))
                    dom_less = dom_less + 1;
                elseif (x(i,n_var + k) == x(j,n_var + k))  
                    dom_equal = dom_equal + 1;
                else
                    dom_more = dom_more + 1;
                end
            end
            if dom_less == 0 && dom_equal ~= n_obj
                individual(i).n = individual(i).n + 1;
            elseif dom_more == 0 && dom_equal ~= n_obj
                individual(i).p = [individual(i).p j];
            end
        end   
        if individual(i).n == 0
            x(i,n_obj + n_var + 1) = 1;
            F(front).f = [F(front).f i];
            num_in_front(front,1)=length(F(front).f);% the number of particles in each front 
        end
    end
% Find the subsequent fronts
    while ~isempty(F(front).f)
       Q = [];
       for i = 1 : length(F(front).f)
           if ~isempty(individual(F(front).f(i)).p)
                for j = 1 : length(individual(F(front).f(i)).p)
                    individual(individual(F(front).f(i)).p(j)).n = ...
                        individual(individual(F(front).f(i)).p(j)).n - 1;
                    if individual(individual(F(front).f(i)).p(j)).n == 0
                        x(individual(F(front).f(i)).p(j),n_obj + n_var + 1) = ...
                            front + 1;
                        Q = [Q individual(F(front).f(i)).p(j)];
                    end
               end
           end
       end
       front =  front + 1;
       F(front).f = Q;
       num_in_front(front,1)=length(F(front).f);% the number of particles in each front 
    end
       % Sort the population according to the front number
        [~,index_of_fronts] = sort(x(:,n_obj + n_var + 1));
        for i = 1 : length(index_of_fronts)
            sorted_based_on_front(i,:) = x(index_of_fronts(i),:);
        end 
    current_index = 0;
    %%  Crowding Distance
    selected_pop=[];
    deleted_pop=[];
    crowd_dist_obj = 0;
    y = sorted_based_on_front;
   % previous_index = current_index + 1;
    % for i = 1 : length(x)
     %       y(i,:) = sorted_based_on_front(current_index + i,:);%put the front_th rank into y
     %end
    %current_index = current_index + i;
    for i = n_var+1 : n_obj+n_var
            [~, index_of_objectives] = sort(y(:,i));
            sorted_based_on_objective = [];
            for j = 1 : length(index_of_objectives)
                sorted_based_on_objective(j,:) = y(index_of_objectives(j),:);
            end
            f_max = ...
                sorted_based_on_objective(length(index_of_objectives), i);
            f_min = sorted_based_on_objective(1,  i);

            if length(index_of_objectives)==1
                y(index_of_objectives(1),n_obj + n_var + 1 + i) = 1;  %If there is only one point in current front

            else
                % deal with boundary points in both decision and objective space
                % twice the distance between the boundary points and its nearest neibohood 
                if (f_max - f_min == 0) % only one point in the current Front
                    y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)=rand;
                    y(index_of_objectives(1),n_obj + n_var + 1 + i)=rand;
                else
                    y(index_of_objectives(length(index_of_objectives)),n_obj + n_var + 1 + i)...
                        = 2*(sorted_based_on_objective(length(index_of_objectives), i)-...
                    sorted_based_on_objective(length(index_of_objectives) -1, i))/(f_max - f_min);
                     y(index_of_objectives(1),n_obj + n_var + 1 + i)=2*(sorted_based_on_objective(2, i)-...
                    sorted_based_on_objective(1, i))/(f_max - f_min);
                end
            end
             for j = 2 : length(index_of_objectives) - 1
                next_obj  = sorted_based_on_objective(j + 1, i);
                previous_obj  = sorted_based_on_objective(j - 1,i);
                if (f_max - f_min == 0) % only one point in the current Front
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = rand;
                else
                    y(index_of_objectives(j),n_obj + n_var + 1 + i) = ...
                         (next_obj - previous_obj)/(f_max - f_min);
                end
             end
    end
     crowd_dist_obj = [];
        crowd_dist_obj(:,1) = zeros(size(y,1),1);
        for i = 1 : n_obj
            crowd_dist_obj(:,1) = crowd_dist_obj(:,1) + y(:,n_obj + n_var + 1+n_var + i);
        end
        

        avg_crowd_dist_obj=mean(crowd_dist_obj);
        y(:,n_obj+n_var+2)= crowd_dist_obj;

            %% calculate the distance in decision sapce
        crowd_dist_var = [];
            newpop=y;
            k1=6;
            for p_i=1:size(y,1)
             poptemp=repmat(y(p_i,1:n_var),size(newpop,1),1);
        %      distance(:,p_i)=sqrt(sum((newpop(:,1:n_var)-poptemp(:,1:n_var)).*(newpop(:,1:n_var)-poptemp(:,1:n_var)),2));%the first column is distance from all particle to the first particle
             distance=sqrt(sum((newpop(:,1:n_var)-poptemp(:,1:n_var)).*(newpop(:,1:n_var)-poptemp(:,1:n_var)),2));% popsize*1
             [sorted_dis,ind_dis]=sort(distance');% sorting the distance to p_i th particle in desend order
             

             if length(distance)>k1
                 crowd_dist_var(p_i,1)=sum([1:k1].*sorted_dis(2:k1+1))/k1;%Eq.(8) in Combining Crowding Estimatioyn in Objective and Decision Space With Multiple Selection and Search Strategies for Multi-Objective Evolutionary Optimization

             else
                 crowd_dist_var(p_i,1)=sum([2:length(distance)].*sorted_dis(2:end))/length(distance);% the

             end
             
             
            end
            %crowd_dist_var=crowd_dist_var./max(crowd_dist_var);
             y(:,n_obj+n_var+3)= crowd_dist_var;
            avg_crowd_dist_var=mean(crowd_dist_var);
            objvar=round(avg_crowd_dist_var/avg_crowd_dist_obj);
           
            special_crowd_dist=zeros(size(crowd_dist_obj,1),1);
            
            for i_SCD = 1 : size(crowd_dist_obj,1)
                        special_crowd_dist(i_SCD)=objvar*(1-c1)*crowd_dist_obj(i_SCD)+ c1*crowd_dist_var(i_SCD); % Eq. (7) in the paper     
            end
            y(:,n_obj+n_var+6)=special_crowd_dist;%put the crowding distace in the decision space in the column (n_obj+n_var+4)
             [~,index_sorted_based_SCD]=sort(y(:,n_obj+n_var+6),'descend');
             y=y(index_sorted_based_SCD,:);
              [~,index_of_fronts] = sort(y(:,n_obj + n_var + 1));
             y=y(index_of_fronts,:);
            population=y;
end 