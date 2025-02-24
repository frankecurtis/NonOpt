function plotQPSolverResults

% Copyright (C) 2025 Frank E. Curtis
%
% This code is published under the MIT License.
%
% Author(s) : Frank E. Curtis, Lara Zebiane

% Load data
data = readtable('times.txt');

% Print total time for entire experiment
fprintf('Total time for experiment (s): %f\n',sum(data.Var6));

% Loop over data to include
for set = -1:5

  % Find unique n values
  n_values = unique(data.Var2);

  % Set m functions
  m_funcs = {@(n) (n+1), @(n) (1.5*n), @(n) (2*n)};

  % Parse data by algorithm
  DAS_data = data(data.Var4 == 0,:);
  IPM_data = data(data.Var4 == 1,:);

  % Parse data further by d_* values (for sets 0 through 2) or m values (for sets 3 through 5)
  if 0 <= set && set <= 2
    DAS_data = DAS_data(DAS_data.Var1 == set,:);
    IPM_data = IPM_data(IPM_data.Var1 == set,:);
  elseif set > 2
    DAS_data_all = DAS_data; DAS_data = [];
    IPM_data_all = IPM_data; IPM_data = [];
    for j = 1:size(DAS_data_all,1)
      for i = 1:length(n_values)
        if DAS_data_all(j,:).Var2 == n_values(i) && DAS_data_all(j,:).Var3 == m_funcs{set-2}(n_values(i))
          DAS_data = [DAS_data; DAS_data_all(j,:)];
        end
      end
    end
    for j = 1:size(IPM_data_all,1)
      for i = 1:length(n_values)
        if IPM_data_all(j,:).Var2 == n_values(i) && IPM_data_all(j,:).Var3 == m_funcs{set-2}(n_values(i))
          IPM_data = [IPM_data; IPM_data_all(j,:)];
        end
      end
    end
  end

  % Reset plot data
  DAS_plot_data = [];
  IPM_plot_data = [];

  % Loop over n values
  for i = 1:length(n_values)

    % Grab data
    DAS_plot_data(:,i) = DAS_data.Var6(DAS_data.Var2 == n_values(i),:);
    IPM_plot_data(:,i) = IPM_data.Var6(IPM_data.Var2 == n_values(i),:);

  end

  % Create figure
  f = figure(set+2);

  % Plot DAS
  subplot(1,2,1);
  boxplot(DAS_plot_data);
  xticklabels(string(n_values));
  xtickangle(45);
  grid on;
  xlabel('n','fontsize',14);
  ylabel('CPU time (seconds)','fontsize',14);
  ylim([0 max(1.1*max(max(DAS_plot_data)),1.1*max(max(IPM_plot_data)))]);
  title('DAS','fontsize',14);

  % Plot IPM
  subplot(1,2,2);
  boxplot(IPM_plot_data);
  xticklabels(string(n_values));
  xtickangle(45);
  grid on;
  xlabel('n','fontsize',14);
  ylabel('CPU time (seconds)','fontsize',14);
  ylim([0 max(1.1*max(max(DAS_plot_data)),1.1*max(max(IPM_plot_data)))]);
  title('IPM','fontsize',14);

  % Set figure rendering size
  f.Position = [100 100 1000 400];

  % Save figure
  saveas(f,sprintf('das_ipm_%d.png',set+1))

end