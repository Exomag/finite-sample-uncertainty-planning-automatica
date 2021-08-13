function generatePlotsCaseStudy2()
    % Parameters
    params = loadParams();

    % Data
    N_trials = 10;
    Nsamples = 1e3 : 1e3 : 3e4;
    methods = {'ExactMoments', 'MomentsRobust', 'Scenario'};
    solutions = cell(3, length(Nsamples));

    % Trials
    for indSample = 1 : length(Nsamples)
        params.Ns = Nsamples(indSample);
        fprintf('Trial %d/%d\n', indSample, length(Nsamples));
        for trial = 1 : N_trials
            % Initialize
            rng(trial);
            tempProblemData = getProblemData(params);
            problemData(trial) = tempProblemData;

            % Exact Moments Approach
            method = methods{1};
            firstSolution = solveProblem(params, problemData(trial), method);
            tempSolution = allocateRisk(params, problemData(trial), method, firstSolution);
            solutions{1,indSample}(trial) = tempSolution;

            % Moments Robust Approach
            method = methods{2};
            firstSolution = solveProblem(params, problemData(trial), method);
            tempSolution = allocateRisk(params, problemData(trial), method, firstSolution);
            solutions{2,indSample}(trial) = tempSolution;

            % Scenario Approach
            method = methods{3};
            tempSolution = solveProblem(params, problemData(trial), method);
            solutions{3,indSample}(trial) = tempSolution;
        end
    end
    
    % Figure - Costs
    avgCosts = zeros(3, length(Nsamples));
    for i = 1 : 3
        for j = 1 : length(Nsamples)
            for k = 1 : N_trials
                avgCosts(i,j) = avgCosts(i,j) + solutions{i,j}(k).objVal;
            end
            avgCosts(i,j) = avgCosts(i,j) / N_trials;
        end
    end
    h1 = figure(); hold on; grid on;
    plot(Nsamples, avgCosts(2,:), 'b-o');
    semilogx(Nsamples, avgCosts(1,:), '-s', 'Color', [0 100/255 0]);
    plot(Nsamples, avgCosts(3,:), 'm-d');
    xlim([Nsamples(1) Nsamples(end)]); xlabel('number of samples $$N_s$$', 'interpreter', 'latex');
    ylabel('total cost', 'interpreter',' latex');
    legend({'MRA', 'EMA', 'SA'}, 'Location', 'NorthWest', 'interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    MagInset(h1, -1, [2e4 3e4 1.8 2], [1.75e4 2.75e4 2.2 2.9], {'NW','SW';'NE','SE'}); legend('hide'); set(gca, 'ytick', get(gca, 'ylim'));
    save2tikz('plots/OpenLoop_Costs');

    % Figure - Solve times
    avgSolveTimes = zeros(3, length(Nsamples));
    for i = 1 : 3
        for j = 1 : length(Nsamples)
            for k = 1 : N_trials
                avgSolveTimes(i,j) = avgSolveTimes(i,j) + solutions{i,j}(k).diagnostic.solvertime;
            end
            avgSolveTimes(i,j) = avgSolveTimes(i,j) / N_trials;
        end
    end
    h2 = figure(); hold on; grid on;
    plot(Nsamples, avgSolveTimes(2,:), 'b-o');
    plot(Nsamples, avgSolveTimes(1,:), '-s', 'Color', [0 100/255 0]);
    plot(Nsamples, avgSolveTimes(3,:), 'm-d');
    xlim([Nsamples(1) Nsamples(end)]); xlabel('number of samples $$N_s$$', 'interpreter', 'latex');
    ylabel('solver time (\si{\second})', 'interpreter',' latex');
    legend({'MRA', 'EMA', 'SA'}, 'Location', 'NorthWest', 'interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    MagInset(h2, -1, [2e4 3e4 0.3 0.6], [1.75e4 2.75e4 3 10], {'NW','SW';'NE','SE'}); legend('hide');  set(gca, 'ytick', get(gca, 'ylim'));
    save2tikz('plots/OpenLoop_SolveTimes');
end
