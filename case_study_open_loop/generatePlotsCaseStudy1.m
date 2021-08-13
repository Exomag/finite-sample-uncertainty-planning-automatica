function generatePlotsCaseStudy1()
    % Parameters
    params = loadParams();

    % Data
    N_trials = 100;
    methods = {'ExactMoments', 'MomentsRobust', 'Scenario'};
    solutions = cell(3, 1);

    % Trials
    for trial = 1 : N_trials
        % Initialize
        rng(trial);
        tempProblemData = getProblemData(params);
        problemData(trial) = tempProblemData;
        fprintf('Trial %d/%d\n', trial, N_trials);

        % Exact Moments Approach
        method = methods{1};
        firstSolution = solveProblem(params, problemData(trial), method);
        tempSolution = allocateRisk(params, problemData(trial), method, firstSolution);
        tempSolution.pViol = runMonteCarlo(params, problemData(trial), tempSolution);
        solutions{1}(trial) = tempSolution;

        % Moments Robust Approach
        method = methods{2};
        firstSolution = solveProblem(params, problemData(trial), method);
        tempSolution = allocateRisk(params, problemData(trial), method, firstSolution);
        tempSolution.pViol = runMonteCarlo(params, problemData(trial), tempSolution);
        solutions{2}(trial) = tempSolution;

        % Scenario Approach
        method = methods{3};
        tempSolution = solveProblem(params, problemData(trial), method);
        tempSolution.pViol = runMonteCarlo(params, problemData(trial), tempSolution);
        solutions{3}(trial) = tempSolution;
    end

    % Figure - Costs
    costsDataBoxplot = [[solutions{2}.objVal]',  [solutions{1}.objVal]', [solutions{3}.objVal]'];
    figure(); hold on;
    boxplot(costsDataBoxplot, 'Labels', {'MRA', 'EMA', 'SA'}, 'Whisker', Inf);
    ylabel('total cost', 'interpreter',' latex');
    legend('hide');
    set(gca, 'TickLabelInterpreter', 'latex');
    save2tikz('plots/OpenLoop_CostBoxplot');

    % Figure - Violations
    violsDataBoxplot = 100 * [[solutions{2}.pViol]', [solutions{1}.pViol]', [solutions{3}.pViol]'];
    epsilon = 100 * params.epsilon;
    figure(); hold on;
    boxplot(violsDataBoxplot, 'Labels', {'MRA', 'EMA', 'SA'}, 'Whisker', Inf);
    xl = xlim;
    plot([xl(1) xl(end)], [epsilon epsilon], 'k--');
    ylim([0 epsilon*1.1]); ylabel('violation probability (\si{\percent})', 'interpreter',' latex');
    legend('hide');
    set(gca, 'TickLabelInterpreter', 'latex');
    save2tikz('plots/OpenLoop_ProbabilityBoxplot');

    %% Figure - Trajectories
    demoInd = 1;
    h3 = figure(); hold on; grid on;
    plot(solutions{2}(demoInd).xVal(1,:), solutions{2}(demoInd).xVal(2,:), 'b-o');
    plot(solutions{1}(demoInd).xVal(1,:), solutions{1}(demoInd).xVal(2,:), '-s', 'Color', [0 100/255 0]);
    plot(solutions{3}(demoInd).xVal(1,:), solutions{3}(demoInd).xVal(2,:), 'm-d');
    plot(problemData(demoInd).xTarget, problemData(demoInd).yTarget, 'rx', 'MarkerSize', 12);
    plot(params.wallPolygonVertMean(:,1), params.wallPolygonVertMean(:,2), 'k', 'LineWidth', 1);
    for j = 1 : 10
        wallPolygonVert = getWallPolygonVert({problemData(demoInd).wallCoefsSamples{1,1}(j,:), ...
            problemData(demoInd).wallCoefsSamples{1,2}(j,:)}, params.yMin, params.xMax);
        plot(wallPolygonVert(:,1), wallPolygonVert(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.25);
    end
    xlim([params.xMin, params.xMax]); xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex');
    ylim([params.yMin, params.yMax]); ylabel('$$x_2$$ (\si{\meter})', 'interpreter', 'latex');
    legend({'MRA', 'EMA', 'SA'}, 'Location', 'NorthWest', 'interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    MagInset(h3, -1, [6.25 6.75 6.8 7.2], [4.5 7.5 1.5 3.5], {'SW','NW';'SE','NE'}); grid off; legend('hide'); ...
        set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    save2tikz('plots/OpenLoop_Trajectories');
end
