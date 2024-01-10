function [val] = GaussQuad(f, a, b, num_points)
  
  if num_points==2
    % 2 point Gauss
    fprintf('using 2 points Gauss Quad\n')
    % weights
    w = ones(2,1); % same as w = [1; 1]
    % gauss points
    xi = zeros(2,1);
    xi(1) = 1/sqrt(3);
    xi(2) = -xi(1);
  elseif num_points==5
    % 5 point Gauss
    fprintf('using 5 points Gauss Quad\n')
    % weights (Based on the provided table in the question)
    w = zeros(5,1);
    w(1) = 0.236926885056189;
    w(2) = 0.478628670499366;
    w(3) = 0.568888888888889;
    w(4) = 0.478628670499366;
    w(5) = 0.236926885056189;
    % gauss points
    xi = zeros(5,1);
    xi(1) = 0.906179845938664;
    xi(2) = 0.538469310105683;
    xi(3) = 0.000000000000000;
    xi(4) = -0.538469310105683;
    xi(5) = -0.906179845938664;
  else
    fprintf('did not enter 2 or 5\n')
  end
    %calculate jacobian
    J=(b-a)/2;
    %map xi to domain of integration
    x=J*xi+(b+a)/2;
    %perform quadrature
    val=J*sum(w .* f(x));
    
  end