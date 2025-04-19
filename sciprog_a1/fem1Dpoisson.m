
% Due to March, 18, 2025
function fem1dpoisson()
% solves -u'' = 1 on (0,1) with u(0) = u(1) = 0 using FEM
% uses piecewise linear basis funcs on a uniform mesh
% plots both the energy norm error and an example solution

    % define the list of element counts for different refinements
    number_of_elements_list = [4, 8, 16, 32, 64, 128];
    energy_norm_errors = zeros(size(number_of_elements_list)); % store energy norm errors
    mesh_widths = 1 ./ number_of_elements_list; % compute mesh widths

    % loop over each refinement level
    for index = 1:length(number_of_elements_list)
        num_elements = number_of_elements_list(index);
        mesh_spacing = 1 / num_elements;

        [stiffness_matrix, load_vector] = assemble_system(num_elements, mesh_spacing);
        fem_interior_solution = stiffness_matrix \ load_vector;

        % include bc: (u(0)=0 and u(1)=0)
        complete_solution = [0; fem_interior_solution; 0];

        % compute the energy norm error for this mesh
        energy_norm_errors(index) = compute_energy_error(complete_solution, mesh_spacing);
    end


%% 


%_________________________________plot energy norm error versus mesh width on a log-log scale
    figure;
    loglog(mesh_widths, energy_norm_errors, 'o-', 'lineWidth', 1, 'markerSize', 8);
    grid on;
    xlabel('mesh width, h');
    ylabel('energy norm error');
    title('energy norm error vs. mesh width for 1d poissan');


%_________________________________approx solution versus the exact solution for one mesh
    
    sample_elements = 10;
    sample_spacing = 1 / sample_elements;
    [sample_stiffness, sample_load] = assemble_system(sample_elements, sample_spacing);
    sample_interior_solution = sample_stiffness \ sample_load;
    sample_solution = [0; sample_interior_solution; 0];

    % mesh nodes for the sample solution, [0,0.1,0.2..1]
    x_nodes = linspace(0, 1, sample_elements + 1);

    % grid to plot the exact solution
    x_fine = linspace(0, 1, 200);
    u_exact = exact_solution(x_fine);

    figure;
    plot(x_nodes, sample_solution, 'bo-', 'lineWidth', 1, 'markerSize', 6, 'DisplayName', 'fem approx');
    hold on;
    plot(x_fine, u_exact, 'r-', 'lineWidth', 1, 'DisplayName', 'exact solution');
    xlabel('x');
    ylabel('u(x)');
    legend('location', 'best');
    title(sprintf('fem vs. exact solution (elements = %d)', sample_elements));
    grid on;
end


%% 


%_________________________________function to assemble the stiffness matrix and load vector
function [stiffness_matrix, load_vector] = assemble_system(num_elements, mesh_spacing)
% assemble the finite element system on the interval [0,1], considering only the interior nodes
% the stiffness matrix is tridiagonal
%   main diag: 2 / mesh_spacing
%   off diag: -1 / mesh_spacing
% the load vector is computed with f(x) = 1, so each entry equals mesh_spacing

    main_diagonal=2*ones(num_elements - 1, 1);
    off_diagonal=-1*ones(num_elements - 2, 1);

    stiffness_matrix = (1/mesh_spacing) * (diag(main_diagonal) + diag(off_diagonal, 1) + diag(off_diagonal, -1));
    load_vector = mesh_spacing * ones(num_elements - 1, 1);
end




%% 



%_________________________________function to compute the energy norm error


function error_value = compute_energy_error(solution_vector, mesh_spacing) 
% calc energy norm error between exact' and FEM'
% basically sqrt of int of (exact' - fem')^2
% integral’s done with 2-pt Gauss quad

    num_intervals = length(solution_vector) - 1;
    sum_squared_error = 0;

    for element = 1:num_intervals
        x_left = (element - 1) * mesh_spacing;
        x_right = element * mesh_spacing;

        fem_derivative = (solution_vector(element + 1) - solution_vector(element))/mesh_spacing;    %the constant fem derivative on this interval

        %use 2 pt gaussian quadrature
        gauss_points = [-1/sqrt(3), 1/sqrt(3)];
        gauss_weights = [1, 1];

        for q = 1:2 % q is 1 or 2
            %map the gaussian point to the current interval
            xi = gauss_points(q); 
            quad_point = (x_right - x_left) / 2 * xi + (x_right + x_left) / 2;
            weight = gauss_weights(q) * (x_right - x_left) / 2; % weight depednds on the interval

            exact_derivative = (1 - 2 * quad_point) / 2;
            difference_value = exact_derivative - fem_derivative;
            sum_squared_error = sum_squared_error + difference_value^2 * weight;
        end
    end

    error_value = sqrt(sum_squared_error);
end

%% 


%_________________________________function to compute the exact solution
function u_exact = exact_solution(x)
    u_exact = 0.5 * x .* (1 - x);
end
