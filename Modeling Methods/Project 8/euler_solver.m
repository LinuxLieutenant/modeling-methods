function [t, y] = euler_solver(f, dt, tf, y0)

    % Find the number of steps, N
    % Round up in MATLAB using ceil() always rounds up to nearest integer
    N = ceil(tf/dt);
    % Discretize time vector, t, using number of steps, N
    t=0:dt:dt*N;
    % Define y vector using initial condition
    y(1)=y0;
    % Loop through time and perform Euler's method
    for i=1:N
        y(i+1)=y(i)+dt*f(t(i),y(i));
    end
end