function [ varargout ] = solveMyProblem( problem,guess,options )
%solveMyProblem - main file for solving NLPs
%
% Syntax:  [ varargout ] = solveMyProblem( problem,guess,options )
%
% Copyright (C) 2019 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the MIT License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.5 
% 1 Aug 2019
% iclocs@imperial.ac.uk


%------------- BEGIN CODE --------------

% Multi-phase problem
if isfield(options,'mp')
    switch options.mp.meshstrategy
        case{'fixed','hp_flexible'}
            [infoNLP,data,options]=transcribeMultiphaseOCP(problem,guess,options);% Format for NLP solver
            [solution,status,data]=solveSingleNLP_DirectCollocation_MultiPhase(infoNLP,data);% Solve the NLP
            [solution] = runPostSolveTasks(problem,solution,options,data); % Output solutions
            if (strcmp(options.resultRep,'res_min_final_manual') || strcmp(options.resultRep,'res_min_final_default'))
                data.options.resultRep='res_min';
                [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions
            end

            varargout{1}=solution;
            varargout{2}=status;
        case{'mesh_refinement'}
            for i=1:length(options.phaseoptions)
                options.phaseoptions{i}.constraintErrorTol_org=problem.phases{i}.constraintErrorTol;
            end
            
            nphase=length(problem.phases);
            errorHistory=cell(1,nphase);
            ConstraintErrorHistory=cell(1,nphase);
            timeHistory=zeros(1,1);
            iterHistory=zeros(1,1);
            solutionHistory=cell(1,1);
            statusHistory=cell(1,1);
            if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                resErrorHistory=cell(1,nphase);
            end

            
            runCondition=1;
            i=1; imax=options.mp.maxMRiter;

            problem_iter=problem;
            while runCondition

                if i~=1 && isfield(options,'ECH') && options.ECH.enabled
                    [ problem_iter,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, i );
                end

                [infoNLP,data,options]=transcribeMultiphaseOCP(problem_iter,guess,options); % Format for NLP solver

                if isfield(options.mp,'regstrategy') && strcmp(options.mp.regstrategy,'simultaneous')
                    if isfield(data.mpdata.data.penalty,'i') && isfield(data.mpdata.data.penalty,'values')
                        idx_penalty=i;
                        idx_penalty(idx_penalty>length(data.mpdata.data.penalty.values))=length(data.mpdata.data.penalty.values);
                        data.mpdata.data.penalty.i=idx_penalty;
                        for j=1:length(options.phaseoptions)
                            data.phasedata{j}.data.penalty.i=i;
                        end
                    else
                        error('Regularization Parameters Not Properly Configured!')
                    end
                end

                [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
                [solution] = runPostSolveTasks(problem,solution,options,data);    % Output solutions

                for j=1:nphase
                    errorHistory{i,j}=max(abs(solution.phaseSol{j}.Error));
                    ConstraintErrorHistory{i,j}=max(solution.phaseSol{j}.ConstraintError);
                    if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                        resErrorHistory{i,j}=solution.phaseSol{j}.residuals.r;
                    end
                    
                end
                timeHistory(i)=solution.mp.computation_time;
                solutionHistory{i}=solution;
                statusHistory{i}=status;

                runCondition_MR=0;
                switch options.mp.errortype
                case{'local_abs'}
                    for j=1:nphase
                        runCondition_MR= (runCondition_MR || any(errorHistory{i,j}>problem.phases{j}.states.xErrorTol_local) || any(ConstraintErrorHistory{i,j}>problem.phases{j}.constraintErrorTol)) && i<=imax;
                    end
                    for j=1:nphase
                        if ~runCondition_MR && (strcmp(options.phaseoptions{j}.resultRep,'res_min_final_manual') || strcmp(options.phaseoptions{j}.resultRep,'res_min_final_default'))
                            data.phasedata{j}.options.resultRep='res_min';
                            [solution] = runPostSolveTasks(problem,solution,options,data);         % Output solutions
                            errorHistory{i,j}=max(abs(solution.phaseSol{j}.Error));
                            ConstraintErrorHistory{i,j}=max(solution.phaseSol{j}.ConstraintError);
                            if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                                resErrorHistory{i,j}=solution.phaseSol{j}.residuals.r;
                            end
                        end
                        timeHistory(i)=solution.mp.computation_time;
                        solutionHistory{i}=solution;
                        statusHistory{i}=status;
                    end

                end
                runCondition=runCondition_MR;
                if runCondition_MR
                    for j=1:nphase
                        [ options.phaseoptions{j}, guess.phases{j} ] = doMeshRefinement( options.phaseoptions{j}, problem.phases{j}, guess.phases{j}, data.phasedata{j}, solution.phaseSol{j}, i );
                    end
                else
                    if isfield(options.mp,'regstrategy') && strcmp(options.mp.regstrategy,'simultaneous') && data.mpdata.data.penalty.i<length(data.mpdata.data.penalty.values)
                        for j=1:nphase
                            [ options.phaseoptions{j}, guess.phases{j}] = doWarmStart( options.phaseoptions{j}, guess.phases{j}, solution.phaseSol{j}, data.phasedata{j} );
                        end
                        runCondition=1;
                    end
                end
                i=i+1;
            end

            if isfield(data.mpdata.data,'penalty') && strcmp(options.mp.regstrategy,'MR_priority')
                if isfield(data.mpdata.data.penalty,'i') && isfield(data.mpdata.data.penalty,'values')
                    for j=1:length(data.mpdata.data.penalty.values)
                        data.mpdata.data.penalty.i=j;
                        for k=1:length(options.phaseoptions)
                            data.phasedata{k}.data.penalty.i=j;
                        end
                        [solution,status,data]=solveSingleNLP_DirectCollocation_MultiPhase(infoNLP,data);
                        [solution] = runPostSolveTasks(problem,solution,options,data);    % Output solutions

                        for k=1:nphase
                            errorHistory{i,k}=max(abs(solution.phaseSol{k}.Error));
                            ConstraintErrorHistory{i,k}=max(solution.phaseSol{k}.ConstraintError);
                            if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                                resErrorHistory{i,k}=solution.phaseSol{k}.residuals.r;
                            end
                        end
                        timeHistory(i)=solution.mp.computation_time;
                        solutionHistory{i}=solution;
                        statusHistory{i}=status;


                        for k=1:length(solution.phaseSol)
                            data.phasedata{k}.multipliers.lambda=solution.phaseSol{k}.multipliers.lambdaNLP;
                        end
                        infoNLP.mpinfoNLP.z0=solution.mp.z_org;
                        i=i+1;
                    end
                else
                    error('Regularization Parameters Not Properly Configured!')
                end
            end

            MeshRefinementHistory.errorHistory=errorHistory;
            MeshRefinementHistory.timeHistory=timeHistory;
            MeshRefinementHistory.iterHistory=iterHistory;
            MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;
            if isfield(options.mp.print,'residual_error') && options.mp.print.residual_error
                MeshRefinementHistory.resErrorHistory=resErrorHistory;
            end
            MeshRefinementHistory.statusHistory=statusHistory;

            varargout{1}=solution;
            varargout{2}=MeshRefinementHistory;
        otherwise
            error('Unknown Meshing Strategy Selected!')
    end

    
    
    
    
else % single phase problem
    switch options.meshstrategy
        case{'fixed','hp_flexible'}
            
            if strcmp(options.transcription,'integral_res_min') && strcmp(options.min_res_mode,'weightedCost')
                if isfield(options,'resCostWeight') 
                    for i=1:length(options.resCostWeight)
                        if isfield(problem.data,'resNormCusWeight')
                            problem.data.resNormCusWeight=problem.data.resNormCusWeight*options.resCostWeight(i);
                        else
                            problem.data.resNormCusWeight=options.resCostWeight(i);
                        end
                        [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
                        [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
                        [ options, guess] = doWarmStart( options, guess, solution, data.dataNLP );
                    end
                    [solution]=runPostSolveTasks(problem, solution,options,data);          % Output solutions

                    varargout{1}=solution;
                    varargout{2}=status;
                else
                    error('Weighting Parameters Not Properly Configured for Residual Minimization with Penalty Method!')
                end
            else
                [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
                [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
                [solution]=runPostSolveTasks(problem, solution,options,data);          % Output solutions
                if (strcmp(options.resultRep,'res_min_final_manual') || strcmp(options.resultRep,'res_min_final_default'))
                    data.options.resultRep='res_min';
                    [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions
                end

                varargout{1}=solution;
                varargout{2}=status;
            end

            
        case{'mesh_refinement'}
            if strcmp(options.transcription,'integral_res_min')
                error('Mesh refinement does not currectly supported the integrated residual minimization method. Please use a fixed mesh for integrated residual minimization, or use direct collocation method with mesh refinement.')
            end
            options.constraintErrorTol_org=problem.constraintErrorTol;
            errorHistory=cell(1,1);
            ConstraintErrorHistory=cell(1,1);
            timeHistory=zeros(1,1);
            iterHistory=zeros(1,1);
            solutionHistory=cell(1,1);
            problemHistory=cell(1,1);
            statusHistory=cell(1,1);
            resErrorHistory=cell(1,1);
            MRiterCheck=false(1);
            
            runCondition=1;
            i=1; imax=options.maxMRiter;

            problem_iter=problem;
            while runCondition

                if isfield(options,'ECH') && options.ECH.enabled
                    if i~=1
                        [ problem_iter,guess,options ] = selectAppliedConstraint( problem, guess, options, data, solutionHistory, i );
                    end
                    problemHistory{i}=problem_iter;
                end

                [infoNLP,data,options]=transcribeOCP(problem_iter,guess,options); % Format for NLP solver

                if isfield(options,'regstrategy') && strcmp(options.regstrategy,'simultaneous')
                    if isfield(data.data.penalty,'i') && isfield(data.data.penalty,'values')
                        idx_penalty=i;
                        idx_penalty(idx_penalty>length(data.data.penalty.values))=length(data.data.penalty.values);
                        data.data.penalty.i=idx_penalty;
                    else
                        error('Regularization Parameters Not Properly Configured!')
                    end
                end

                [solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
                [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions

                maxAbsError=max(abs(solution.Error));
                maxAbsConstraintError=max(solution.ConstraintError);
                if isfield(options.print,'residual_error') && options.print.residual_error
                    resError=solution.residuals.r;
                    resErrorHistory{i,1}=resError;
                end
                errorHistory{i,1}=maxAbsError;
                ConstraintErrorHistory{i,1}=maxAbsConstraintError;
                timeHistory(i)=solution.computation_time;
                solutionHistory{i,1}=solution;
                statusHistory{i,1}=status;
                iterHistory(i)=status.iter;


                switch options.errortype
                case{'local_abs'}
                    runCondition_MR=(any(maxAbsError>problem.states.xErrorTol_local) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax;
                    if ~runCondition_MR && (strcmp(options.resultRep,'res_min_final_manual') || strcmp(options.resultRep,'res_min_final_default'))
                        data.options.resultRep='res_min';
                        [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions
                        maxAbsError=max(abs(solution.Error));
                        maxAbsConstraintError=max(solution.ConstraintError);
                        errorHistory{i,1}=maxAbsError;
                        ConstraintErrorHistory{i,1}=maxAbsConstraintError;
                        timeHistory(i)=solution.computation_time;
                        solutionHistory{i,1}=solution;
                        statusHistory{i,1}=status;
                        if isfield(options.print,'residual_error') && options.print.residual_error
                            resError=solution.residuals.r;
                            resErrorHistory{i,1}=resError;
                        end
                        iterHistory(i)=status.iter;
                    end
                    if i>1
                          MRiterCheck(i)=(any((min(cell2mat(errorHistory(1:i-1)))-errorHistory{i})./min(cell2mat(errorHistory(1:i-1)))<0) || all(0<(min(cell2mat(errorHistory(1:i-1)))-errorHistory{i})./min(cell2mat(errorHistory(1:i-1)))<0.05)) && (any((min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0) || all(0<(min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0.05));
                    end
                    if runCondition_MR && i>5 && all(MRiterCheck(i-4:i))
                        waitAnswer=1;
                        while waitAnswer
                            keepMR = input('Possible slow convergence or diverging mesh refinement iterations, continue to refine the mesh? (Yes/No) \n', 's');
                            if strcmp(keepMR, 'Yes')
                                waitAnswer=0;
                            elseif strcmp(keepMR, 'No')
                                runCondition_MR=false(1); 
                                waitAnswer=0;
                            else
                                disp('Answer not recognized, please enter again!')
                            end
                        end
                    end

                case{'int_res'}
                    runCondition_MR=(any(resError*0.99>problem.states.xErrorTol_integral'.^2) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax;
                    if ~any(maxAbsError>problem.states.xErrorTol_local) && ~any(maxAbsConstraintError>problem.constraintErrorTol) && runCondition_MR
                        solution.Error=solution.Error./max(solution.Error).*problem.states.xErrorTol_local.*resError'./problem.states.xErrorTol_integral;
                    end
                    
                    if i>1
                          MRiterCheck(i)=(any((min(cell2mat(resErrorHistory(1:i-1)))-resErrorHistory{i})./min(cell2mat(resErrorHistory(1:i-1)))<0) || all(0<(min(cell2mat(resErrorHistory(1:i-1)))-resErrorHistory{i})./min(cell2mat(resErrorHistory(1:i-1)))<0.05)) && (any((min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0) || all(0<(min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0.05));
                    end
                    if runCondition_MR && i>5 && all(MRiterCheck(i-4:i))
                        waitAnswer=1;
                        while waitAnswer
                            keepMR = input('Possible slow convergence or diverging mesh refinement iterations, continue to refine the mesh? (Yes/No) \n', 's');
                            if strcmp(keepMR, 'Yes')
                                waitAnswer=0;
                            elseif strcmp(keepMR, 'No')
                                runCondition_MR=false(1); 
                                waitAnswer=0;
                            else
                                disp('Answer not recognized, please enter again!')
                            end
                        end
                    end
                    
                case{'both'}
                    runCondition_local=(any(maxAbsError>problem.states.xErrorTol_local) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax;
                    runCondition_integral=(any(resError*0.99>problem.states.xErrorTol_integral'.^2) || any(maxAbsError>problem.states.xErrorTol_local) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax;
                    runCondition_MR=runCondition_local || runCondition_integral;
                    if runCondition_integral && ~runCondition_local && (strcmp(options.resultRep,'res_min_final_manual') || strcmp(options.resultRep,'res_min_final_default'))
                        data.options.resultRep='res_min';
                        [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions
                        maxAbsError=max(abs(solution.Error));
                        maxAbsConstraintError=max(solution.ConstraintError);
                        errorHistory{i,1}=maxAbsError;
                        ConstraintErrorHistory{i,1}=maxAbsConstraintError;
                        timeHistory(i)=solution.computation_time;
                        solutionHistory{i,1}=solution;
                        statusHistory{i,1}=status;
                        if isfield(options.print,'residual_error') && options.print.residual_error
                            resError=solution.residuals.r;
                            resErrorHistory{i,1}=resError;
                        end
                        iterHistory(i)=status.iter;

                       runCondition_MR=(any(resError*0.99>problem.states.xErrorTol_integral'.^2) || any(maxAbsError>problem.states.xErrorTol_local) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax;
                    end
                    
                    if i>1
                          MRiterCheck(i)=((any((min(cell2mat(errorHistory(1:i-1)))-errorHistory{i})./min(cell2mat(errorHistory(1:i-1)))<0) || all(0<(min(cell2mat(errorHistory(1:i-1)))-errorHistory{i})./min(cell2mat(errorHistory(1:i-1)))<0.05)) && (any((min(cell2mat(resErrorHistory(1:i-1)))-resErrorHistory{i})./min(cell2mat(resErrorHistory(1:i-1)))<0) || all(0<(min(cell2mat(resErrorHistory(1:i-1)))-resErrorHistory{i})./min(cell2mat(resErrorHistory(1:i-1)))<0.05))) && (any((min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0) || all(0<(min(cell2mat(ConstraintErrorHistory(1:i-1)))-ConstraintErrorHistory{i})./min(cell2mat(ConstraintErrorHistory(1:i-1)))<0.05));
                    end
                    if runCondition_MR && i>5 && all(MRiterCheck(i-4:i))
                        waitAnswer=1;
                        while waitAnswer
                            keepMR = input('Possible slow convergence or diverging mesh refinement iterations, continue to refine the mesh? (Yes/No) \n', 's');
                            if strcmp(keepMR, 'Yes')
                                waitAnswer=0;
                            elseif strcmp(keepMR, 'No')
                                runCondition_MR=false(1); 
                                waitAnswer=0;
                            else
                                disp('Answer not recognized, please enter again!')
                            end
                        end
                    end
                end
                
                runCondition=runCondition_MR;
                if runCondition_MR
                    [ options, guess ] = doMeshRefinement( options, problem, guess, data, solution, i );
                else
                    if isfield(options,'regstrategy') && strcmp(options.regstrategy,'simultaneous') && data.data.penalty.i<length(data.data.penalty.values)
                        [ options, guess] = doWarmStart( options, guess, solution, data );
                        runCondition=1;
                    end
                end
                i=i+1;


            end

            if isfield(data.data,'penalty') && strcmp(data.options.regstrategy,'MR_priority')
                if isfield(data.data.penalty,'i') && isfield(data.data.penalty,'values')
                    for j=1:length(data.data.penalty.values)
                        data.data.penalty.i=j;
                        [solution,status,data]=solveSingleNLP_DirectCollocation(infoNLP,data);
                        [solution]=runPostSolveTasks(problem,solution,options,data);         % Output solutions

                        maxAbsError=max(abs(solution.Error));
                        maxAbsConstraintError=max(solution.ConstraintError);
                        errorHistory{i,1}=maxAbsError;
                        ConstraintErrorHistory{i,1}=maxAbsConstraintError;
                        timeHistory(i)=solution.computation_time;
                        solutionHistory{i,1}=solution;
                        statusHistory{i,1}=status;
                        if isfield(options.print,'residual_error') && options.print.residual_error
                            resError=solution.residuals.r;
                            resErrorHistory{i,1}=resError;
                        end
                        iterHistory(i)=status.iter;

                        data.multipliers.lambda=solution.multipliers.lambdaNLP;
                        infoNLP.z0=solution.z;
                        i=i+1;
                    end
                else
                    error('Regularization Parameters Not Properly Configured!')
                end
            end

            MeshRefinementHistory.errorHistory=errorHistory;
            MeshRefinementHistory.timeHistory=timeHistory;
            MeshRefinementHistory.iterHistory=iterHistory;
            MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;
            if isfield(options.print,'residual_error') && options.print.residual_error
                MeshRefinementHistory.resErrorHistory=resErrorHistory;
            end
            MeshRefinementHistory.statusHistory=statusHistory;
            MeshRefinementHistory.solutionHistory=solutionHistory;
            if isfield(options,'ECH') && options.ECH.enabled
                MeshRefinementHistory.problemHistory=problemHistory;
            end
            varargout{1}=solution;
            varargout{2}=MeshRefinementHistory;
        otherwise
            error('Unknown Meshing Strategy Selected!')
    end
end
