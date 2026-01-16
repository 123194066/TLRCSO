function [Best_Score,Convergence_curve,BestFit, f1_zong, f2_zong, f3_zong, f4_zong]=BSA(SearchAgents_no,dim,Max_iteration,ub,lb,pos_tag,all_param)

   Tma=Max_iteration-1;   
    FQ = 2;   
    c1 = 2;
    c2 = 2;
    a1 = 1;
    a2 = 1;
%end

% set the parameters
xmi=lb'; % Lower bounds
xma=ub'; % Upper bounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
TOP= initialization(SearchAgents_no, dim);

x=TOP'; 
for i = 1 : SearchAgents_no
    %x( i, : ) = xmi + (xma - xmi) .* rand( 1, dimension ); 
    PP=reshape(x(i,:),2,dim/2)';
     [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
    fit(i) = W;
    f1_x(i) = f1;
    f2_x(i) = f2;
    f3_x(i) = f3;
    f4_x(i) = f4;

end
pFit = fit; % The individual's best fitness value
pX = x;     % The individual's best position corresponding to the pFit

f1Fit = f1_x;
f2Fit = f2_x;
f3Fit = f3_x;
f4Fit = f4_x;
[ fMin, bestIndex ] = min( fit );  % fMin denotes the global optimum
% bestX denotes the position corresponding to fMin
bestX = x( bestIndex, : );   
bestX = x(bestIndex, :);
% 初始化收敛曲线
Convergence_curve = zeros(1, Max_iteration);
f1_zong = zeros(1, Max_iteration);
f2_zong = zeros(1, Max_iteration);
f3_zong = zeros(1, Max_iteration);
f4_zong = zeros(1, Max_iteration);
Convergence_curve(1) = fMin;
f1_zong(1) = f1Fit(bestIndex);
f2_zong(1) = f2Fit(bestIndex);
f3_zong(1) = f3Fit(bestIndex);
f4_zong(1) = f4Fit(bestIndex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the iteration.

 for t = 1 : Tma

     disp(['BSA代数','=',num2str(t)]);
    prob = rand( SearchAgents_no, 1 ) .* 0.2 + 0.8;%The probability of foraging for food

    if( mod( t, FQ ) ~= 0 )         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Birds forage for food or keep vigilance
        
        sumPfit = sum( pFit );
        meanP = mean( pX );
        for i = 1 : SearchAgents_no
            if rand < prob(i)
                x( i, : ) = x( i, : ) + c1 * rand.*(bestX - x( i, : ))+ ...
                    c2 * rand.*( pX(i,:) - x( i, : ) );
            else
                person = randiTabu( 1, SearchAgents_no, i, 1 );


              N2=exp( -(pFit(person) - pFit(i))/(abs( pFit(person)-pFit(i) )...
                        + realmin) * pFit(person)/(sumPfit + realmin) * SearchAgents_no ); 
               if N2>=2
                   N2=2;
               end

             x( i, : ) = x( i, : ) + rand.*(meanP - x( i, : )) * a1 * ...
                    exp( -pFit(i)/( sumPfit + realmin) * SearchAgents_no ) + a2 * ...
                    ( rand*2 - 1) .* ( pX(person,:) - x( i, : ) ) * N2;

            end

            x( i, : ) = Bounds( x( i, : ), xmi, xma );  
            PP=reshape(x(i,:),2,dim/2)';
                [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else
        FL = rand( SearchAgents_no, 1 ) .* 0.4 + 0.5;    %The followed coefficient

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide the bird swarm into two parts: producers and scroungers.
        [ans, minIndex ] = min( pFit );
        [ans, maxIndex ] = max( pFit );
        choose = 0;
        if ( minIndex < 0.5*SearchAgents_no && maxIndex < 0.5*SearchAgents_no )
            choose = 1;
        end
        if ( minIndex > 0.5*SearchAgents_no && maxIndex < 0.5*SearchAgents_no )
            choose = 2;
        end
        if ( minIndex < 0.5*SearchAgents_no && maxIndex > 0.5*SearchAgents_no )
            choose = 3;
        end
        if ( minIndex > 0.5*SearchAgents_no && maxIndex > 0.5*SearchAgents_no )
            choose = 4;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if choose < 3
            for i = (SearchAgents_no/2+1) : SearchAgents_no
                x( i, : ) = x( i, : ) * ( 1 + randn );
                x( i, : ) = Bounds( x( i, : ), xmi, xma );
                    PP=reshape(x(i,:),2,dim/2)';
                [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
            end
            if choose == 1 
                x( minIndex,: ) = x( minIndex,: ) * ( 1 + randn );
                x( minIndex, : ) = Bounds( x( minIndex, : ), xmi, xma );
                    PP=reshape(x(i,:),2,dim/2)';
                 [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
         fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
            end
            for i = 1 : 0.5*SearchAgents_no
                if choose == 2 || minIndex ~= i
                    person = randi( [(0.5*SearchAgents_no+1), SearchAgents_no ], 1 );
                    x( i, : ) = x( i, : ) + (pX(person, :) - x( i, : )) * FL( i );
                    x( i, : ) = Bounds( x( i, : ), xmi, xma );
                    PP=reshape(x(i,:),2,dim/2)';
                [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
         fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
                end
            end
        else
            for i = 1 : 0.5*SearchAgents_no
                x( i, : ) = x( i, : ) * ( 1 + randn );
                x( i, : ) = Bounds( x( i, : ), xmi, xma );
                   PP=reshape(x(i,:),2,dim/2)';
                    [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
        fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
            end
            if choose == 4 
                x( minIndex,: ) = x( minIndex,: ) * ( 1 + randn );
                x( minIndex, : ) = Bounds( x( minIndex, : ), xmi, xma );
                   PP=reshape(x(i,:),2,dim/2)';
                   [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
         fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
            end
            for i = (0.5*SearchAgents_no+1) : SearchAgents_no
                if choose == 3 || minIndex ~= i
                    person = randi( [1, 0.5*SearchAgents_no], 1 );
                    x( i, : ) = x( i, : ) + (pX(person, :) - x( i, : )) * FL( i );
                    x( i, : ) = Bounds( x( i, : ), xmi, xma );
                        PP=reshape(x(i,:),2,dim/2)';
                        [W, f1, f2, f3, f4, ~] = Link_Table(PP, pos_tag, all_param);
         fit(i) = W;
        f1_x(i) = f1;
        f2_x(i) = f2;
        f3_x(i) = f3;
        f4_x(i) = f4;
                end
            end   
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the individual's best fitness vlaue and the global best one

    for i = 1 : SearchAgents_no 
        if fit(i) < pFit(i)
            pFit(i) = fit(i);
            pX(i, :) = x(i, :);
            f1Fit(i) = f1_x(i);
            f2Fit(i) = f2_x(i);
            f3Fit(i) = f3_x(i);
            f4Fit(i) = f4_x(i);
        end

        if( pFit( i ) < fMin )
            fMin = pFit( i );
            bestX = pX( i, : );   
        end
    end
  
     Convergence_curve(t + 1) = fMin;
    [~, bestIndex] = find(pFit == fMin);
    f1_zong(t + 1) = f1Fit(min(bestIndex));
    f2_zong(t + 1) = f2Fit(min(bestIndex));
    f3_zong(t + 1) = f3Fit(min(bestIndex));
    f4_zong(t + 1) = f4Fit(min(bestIndex));

end
Best_Score=fMin;
Convergence_curve=Convergence_curve;
BestFit=bestX';

% End of the main program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions are associated with the main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
