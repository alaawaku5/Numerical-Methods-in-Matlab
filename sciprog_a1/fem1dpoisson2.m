function fem1dpoisson()
% solves -u'' = f on (0,1) with u(0) = u(1) = 0
% targeting a singular solution: u(x) = x^(3/4) * (1 - x)
% mesh isn’t uniform anymore, we’re using x_i = (i/n)^beta

    number_of_elements_list = [8, 16, 32, 64, 128];

    beta = 1.0;

    energy_norm_errors = zeros(size(number_of_elements_list));

    for idx = 1:length(number_of_elements_list)
        num_elems = number_of_elements_list(idx);
       
        xcoords = generate_nonuniform_mesh(num_elems, beta); % non uniform according to beta
        [stiffness_matrix, load_vector] = assemble_system_sing(xcoords); 

        fem_interior_solution = stiffness_matrix \ load_vector;
        
        complete_solution = [0; fem_interior_solution; 0];

        energy_norm_errors(idx) = compute_energy_error_sing(complete_solution, xcoords);
    end


%% 

%_________________________________plot energy norm error versus number of elements list on a log-log scale
    figure;
    loglog(number_of_elements_list, energy_norm_errors, 'o-', 'lineWidth', 1.5, 'markerSize', 8);
    grid on;
    xlabel('number of elements, N');
    ylabel('energy norm error');
    title(sprintf('singular solution: error vs. N (beta = %.2f)', beta));


%_________________________________approx solution versus the exact solution for one mesh
    sample_elems = 16;
    sample_xcoords = generate_nonuniform_mesh(sample_elems, beta);
    [sample_stiff, sample_load] = assemble_system_sing(sample_xcoords);
    sample_interior_sol = sample_stiff \ sample_load;
    sample_full_sol = [0; sample_interior_sol; 0];

    x_nodes = sample_xcoords; % x axis vals

    x_fine = linspace(0,1,400); 
    u_exact_fine = exact_solution_sing(x_fine); % exxact sol

    figure;
    plot(x_nodes, sample_full_sol, 'bo-', 'lineWidth', 1.2, ...
         'markerSize', 6, 'displayName', 'fem approx');
    hold on;
    plot(x_fine, u_exact_fine, 'r-', 'lineWidth', 1.2, ...
         'displayName', 'exact solution');
    xlabel('x');
    ylabel('u(x)');
    legend('location','best');
    title(sprintf('fem approx vs. exact (N=%d, beta=%.2f)', sample_elems, beta));
    grid on;
end


function xcoords = generate_nonuniform_mesh(num_elems, beta)
% xcoords is a vector of length (num_elems+1) with xcoords(1)=0, xcoords(end)=1, and xcoords(i) = (i/num_elems)^beta.
    xcoords = zeros(num_elems+1, 1);
    for i = 0:num_elems
        xcoords(i+1) = (i/num_elems)^beta;
    end
end



%% 

function [stiffness_matrix, load_vector] = assemble_system_sing(xcoords)
% xcoords has length n+1, so there are n = num_elems intervals

    n = length(xcoords) - 1;            % number of elements
    stiffness_matrix = zeros(n-1, n-1); % for interior nodes
    load_vector      = zeros(n-1, 1);

    for elem = 1:n

        left_node = elem; % fix which node we are taking
        right_node = elem + 1;
        x_left = xcoords(left_node);
        x_right = xcoords(right_node);


        h_e = x_right - x_left;

        localK = (1 / h_e) * [1, -1; -1, 1]; % local stiffness (*1)

        x_mid = 0.5*(x_left + x_right); % midpoint (*2)
        localb = h_e * f_sing(x_mid) * [0.5; 0.5];


%____________shifting :  interior node? shift index (ignore boundaries)
        if (left_node >= 2) && (left_node <= n)
            iL = left_node - 1;
        else
            iL = 0;  % boundary
        end


        if (right_node >= 2) && (right_node <= n)
            iR = right_node - 1;
        else
            iR = 0;  % boundary
        end

%____________assemble stiffness matrix 
        % (1,1) -> left node, (2,2) -> right node
        if iL > 0
            stiffness_matrix(iL,iL) = stiffness_matrix(iL,iL) + localK(1,1);
        end
        if iR > 0
            stiffness_matrix(iR,iR) = stiffness_matrix(iR,iR) + localK(2,2);
        end
        if iL>0 && iR>0
            stiffness_matrix(iL,iR) = stiffness_matrix(iL,iR) + localK(1,2);
            stiffness_matrix(iR,iL) = stiffness_matrix(iR,iL) + localK(2,1);
        end

%____________assemble load vector
        if iL > 0
            load_vector(iL) = load_vector(iL) + localb(1);
        end
        if iR > 0
            load_vector(iR) = load_vector(iR) + localb(2);
        end
    end
end




%% 
function val = f_sing(x)
    val = (3/16)*x.^(-5/4).*(1 - x) + (3/2)*x.^(-1/4);
end


%% 

function error_value = compute_energy_error_sing(complete_sol, xcoords)
    n = length(xcoords) - 1;
    sum_sq_error = 0;


    for elem = 1:n
        x_left = xcoords(elem);
        x_right = xcoords(elem+1);
        h_e = x_right - x_left;

        % FEM derivative on [x_left, x_right]
        % complete_sol has n+1 entries (global nodes 1 to n+1)

        fem_deriv = (complete_sol(elem+1) - complete_sol(elem)) / h_e;

        % 2-pt Gauss quadrature to get error in derivatives
        gauss_pts = [-1/sqrt(3), 1/sqrt(3)];
        gauss_wts = [1, 1];

        for q = 1:2
            xi = gauss_pts(q);
            wq = gauss_wts(q);

            % map Gauss pt from [-1,1] to [x_left, x_right]
            xq = 0.5*(x_right - x_left)*xi + 0.5*(x_right + x_left);
            dq = 0.5*(x_right - x_left);

            diff_val = exact_derivative_sing(xq) - fem_deriv;

            sum_sq_error = sum_sq_error + (diff_val^2)*wq*dq;
        end
    end

    error_value = sqrt(sum_sq_error);
end



%% 


function val = exact_solution_sing(x)
    val = x.^(3/4).*(1 - x);
end


function val = exact_derivative_sing(x)
    % u'(x) = (3/4)* x^(-1/4)*(1-x) - x^(3/4).
    val = (3/4)*x.^(-1/4).*(1 - x) - x.^(3/4);
end











