%% Add path

addpath(genpath('MM_testfunctions/functions'));
addpath(genpath('Indicator_calculation/'));
addpath(genpath('MM_testfunctions/'));
addpath(genpath('Indicator_calculation/'));

 clear
  global fname
  N_function=27;% number of test function
  runtimes=30;  % odd number
  fun_num=27;%function number
 %% Initialize the parameters in MMO test functions
     %[2 3]
     for i_func=1: fun_num
        switch i_func
            case 1
                fname='MMF1';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[1 -1];     % the low bounds of the decision variables
                xu=[3 1];      % the up bounds of the decision variables
                repoint=[1.1,1.1]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)     
            case 2
                fname='MMF2';
                n_obj=2;
                n_var=2;
                xl=[0 0];
                xu=[1 2];
                repoint=[1.1,1.1];
            case 3
                fname='MMF4';
                n_obj=2;
                n_var=2;
                xl=[-1 0];
                xu=[1 2];
                repoint=[1.1,1.1];
            case 4
                fname='MMF5';
                n_obj=2;
                n_var=2;
                xl=[1 -1];
                xu=[3 3];
                repoint=[1.1,1.1];
             case 5
                fname='MMF6';
                n_obj=2;
                n_var=2;
                xl=[1 -1];
                xu=[3 2];
                repoint=[1.1,1.1];
            case 6
                fname='MMF7';
                n_obj=2;
                n_var=2;
                xl=[1 -1];
                xu=[3 1];
                repoint=[1.1,1.1];
             case 7
                fname='MMF8';
                n_obj=2;
                n_var=2;
                xl=[-pi 0];
                xu=[pi 9];
               repoint=[1.1,1.1];
              case 8
                fname='MMF9';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.21,11]; % reference point used to calculate the Hypervolume
            case 9
               fname='MMF10';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1];      % the up bounds of the decision variables
               repoint=[1.21,13.2]; % reference point used to calculate the Hypervolume
           case 10
                fname='MMF10_l';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.21,13.2]; % reference point used to cal
          case 11
                fname='MMF11';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.21,15.4];
           case 12
                fname='MMF11_l';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.21,15.4];
   
          case 13
                fname='MMF12';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0 0];     % the low bounds of the decision variables
                xu=[1 1];      % the up bounds of the decision variables
                repoint=[1.54,1.1];
          case 14
                fname='MMF12_l';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[0 0];     % the low bounds of the decision variables
                xu=[1 1];      % the up bounds of the decision variables
                repoint=[1.54,1.1];
         case 15
                 %*need to be modified
                fname='MMF13';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=3;       % the dimensions of the objective space
                xl=[0.1 0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.54,15.4];
         case 16
                 %*need to be modified
                 fname='MMF13_l';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=3;       % the dimensions of the objective space
                xl=[0.1 0.1 0.1];     % the low bounds of the decision variables
                xu=[1.1 1.1 1.1];      % the up bounds of the decision variables
                repoint=[1.54,15.4];
         case 17
                fname='MMF14';  % function name
                n_obj=3;       % the dimensions of the decision space
                n_var=3;       % the dimensions of the objective space
                xl=[0 0 0];     % the low bounds of the decision variables
                xu=[1 1 1];      % the up bounds of the decision variables
                repoint=[2.2,2.2,2.2];
          case 18
                fname='MMF15';  % function name
                n_obj=3;       % the dimensions of the decision space
                n_var=3;       % the dimensions of the objective space
                xl=[0 0 0];     % the low bounds of the decision variables
                xu=[1 1 1];      % the up bounds of the decision variables
                repoint=[2.5,2.5,2.5];
         case 19
                fname='MMF1_z';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[1 -1];     % the low bounds of the decision variables
                xu=[3 1];      % the up bounds of the decision variables
                repoint=[1.1,1.1];
        case 20
                fname='MMF1_e';  % function name
                n_obj=2;       % the dimensions of the decision space
                n_var=2;       % the dimensions of the objective space
                xl=[1 -20];     % the low bounds of the decision variables
                xu=[3 20];      % the up bounds of the decision variables
                repoint=[1.1,1.1];
       case 21
                fname='MMF14_a';  % function name
                n_obj=3;
                n_var=3;
                xl=[0 0 0];
                xu=[1 1 1];
                repoint=[2.2,2.2,2.2];
       case 22
                fname='MMF15_a';  % function name
                n_obj=3;
                n_var=3;
                xl=[0 0 0];
                xu=[1 1 1]; 
                repoint=[2.5,2.5,2.5];
       case 23
                 fname='MMF15_a_l';  % function name
                 n_obj=3;
                 n_var=3;
                 xl=[0 0 0];
                xu=[1 1 1];
                 repoint=[2.5,2.5,2.5];
        case 24 
                   fname='MMF15_l';  % function name
                   n_obj=3;       % the dimensions of the decision space
                    n_var=3;       % the dimensions of the objective space
                    xl=[0 0 0];     % the low bounds of the decision variables
                    xu=[1 1 1];      % the up bounds of the decision variables
                    repoint=[2.5,2.5,2.5];
        case 25
                  fname='MMF16_l1';  % function name
                   n_obj=3;
                   n_var=3;
                    xl=[0 0 0];
                    xu=[1 1 1];
                    repoint=[2.5,2.5,2.5];
   
        case 26
                  fname='MMF16_l2';  % function name
                  n_obj=3;
                  n_var=3;
                  xl=[0 0 0];
                  xu=[1 1 1];
                   repoint=[2.5,2.5,2.5];
        case 27
                fname='MMF16_l3';  % function name
                n_obj=3;
                n_var=3;
                xl=[0 0 0];
                xu=[1 1 1];
                repoint=[2.5,2.5,2.5];
        end
   %% Load reference PS and PF data
         load  (strcat([fname,'_Reference_PSPF_data']));
       %% Initialize the population size and the maximum evaluations
          popsize=100*n_var;
          Max_fevs=5000*n_var;
          Max_Gen=fix(Max_fevs/popsize);
           for j=1:runtimes
               %% Search the PSs using MMODE
               fprintf('Running test function: %s \n %d times \n', fname,j);
               
                  [ps,pf]=SSMMCOASC(fname,xl,xu,n_obj,popsize,Max_Gen,4,1,0.9);
                
            hyp=Hypervolume_calculation(pf,repoint);
              IGDx=IGD_calculation(ps,PS);
              IGDf=IGD_calculation(pf,PF);
              CR=CR_calculation(ps,PS,fname);
              PSP=CR/IGDx ;%;
              hv(i_func,j)=1/hyp;
             psp(i_func,j)=1/PSP; 
            igdx(i_func,j)=IGDx;
            igdf(i_func,j)=IGDf;  
           end
     end 
 
  