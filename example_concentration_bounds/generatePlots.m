function generatePlots()
    % Concentration bounds as a function of confidence level beta
    Ns = 1000;
    betaAx = logspace(-6, 0, 1e4);

    Tval = (3*(Ns-1)) / (Ns-3) * finv(1-betaAx, 3, Ns-3);
    r1_beta = real(sqrt(Tval / Ns));
    r2_beta = max(abs(1 - (Ns-1) ./ chi2inv(betaAx/2,   Ns-1)), ...
                  abs(1 - (Ns-1) ./ chi2inv(1-betaAx/2, Ns-1)));
    figure();
    semilogx(betaAx, r1_beta, '-', 'LineWidth', 1); hold on; grid on;
    semilogx(betaAx, r2_beta, '--', 'LineWidth', 1);
    xlabel('confidence level $$\beta$$', 'interpreter', 'latex');
    ylabel('bound value', 'interpreter', 'latex');
    legend({'$$r_1$$', '$$r_2$$'}, 'Location', 'NorthEast', 'interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    save2tikz('plots/Example_betaDependence');

    % Concentration bounds as a function of number of samples
    NAx = logspace(1, 4, 1e4);
    beta = 1e-3;
    Tval = (3*(NAx-1)) ./ (NAx-3) .* finv(1-beta, 3, NAx-3);
    r1_N = sqrt(Tval ./ NAx);
    r2_N = max(abs(1 - (NAx-1) ./ chi2inv(beta/2,   NAx-1)), ...
               abs(1 - (NAx-1) ./ chi2inv(1-beta/2, NAx-1)));
    figure();
    semilogx(NAx, r1_N, '-', 'LineWidth', 1); hold on; grid on;
    semilogx(NAx, r2_N, '--', 'LineWidth', 1);
    xlabel('number of samples $$N_s$$', 'interpreter', 'latex');
    ylabel('bound value', 'interpreter', 'latex');
    legend({'$$r_1$$', '$$r_2$$'}, 'Location', 'NorthEast', 'interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex');
    save2tikz('plots/Example_sampleDependence');
end
