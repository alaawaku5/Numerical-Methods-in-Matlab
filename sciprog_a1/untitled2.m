import numpy as np
import matplotlib.pyplot as plt

def euler_explicit(f, t0, T, y0, N):
    """
    Implements the explicit (forward) Euler method to solve:
        y'(t) = f(t, y(t)),  y(t0) = y0
    over the interval [t0, T] with N uniform steps.

    Parameters
    ----------
    f   : function
          Right-hand side of the ODE, f(t, y) -> array-like
    t0  : float
          Initial time
    T   : float
          Final time
    y0  : array-like
          Initial state (vector)
    N   : int
          Number of steps

    Returns
    -------
    t : ndarray of shape (N+1,)
        Time points
    y : ndarray of shape (N+1, len(y0))
        Approximate solution at each time step
    """
    # Step size
    h = (T - t0) / N
    
    # Initialize arrays to store time and solution
    t = np.zeros(N+1)
    y = np.zeros((N+1, len(y0)))
    
    # Initial conditions
    t[0] = t0
    y[0] = y0
    
    # Perform Euler steps
    for k in range(N):
        t[k+1] = t[k] + h
        y[k+1] = y[k] + h * f(t[k], y[k])
    
    return t, y

def predator_prey(t, y, r=0.5, K=50.0, alpha=0.02, beta=0.5, delta=0.4):
    """
    Predator-prey system with logistic growth for the prey.

    y = [y1, y2]
        y1' = r * y1 * (1 - y1/K) - alpha * y1 * y2
        y2' = beta * alpha * y1 * y2 - delta * y2

    Parameters
    ----------
    t     : float
            Current time (not used in these equations, but included for generality)
    y     : array-like, shape (2,)
            [y1, y2] at current time
    r     : float
            Intrinsic growth rate of prey
    K     : float
            Carrying capacity of prey
    alpha : float
            Predation rate coefficient
    beta  : float
            Efficiency of converting eaten prey into predator growth
    delta : float
            Natural death rate of predators

    Returns
    -------
    dydt : ndarray of shape (2,)
           Time derivatives [dy1/dt, dy2/dt]
    """
    y1, y2 = y
    dy1 = r * y1 * (1 - y1 / K) - alpha * y1 * y2
    dy2 = beta * alpha * y1 * y2 - delta * y2
    return np.array([dy1, dy2])

# Example usage
if __name__ == "__main__":
    # Parameters
    r = 0.5
    K = 50.0
    alpha = 0.02
    beta = 0.5
    delta = 0.4

    # Time setup
    t0 = 0.0   # start time
    T = 200.0  # end time
    N = 2000   # number of Euler steps

    # Initial populations
    y0 = np.array([10.0, 5.0])  # [prey, predator]

    # Solve using Euler's method
    t_vals, y_vals = euler_explicit(
        lambda t, y: predator_prey(t, y, r, K, alpha, beta, delta),
        t0, T, y0, N
    )

    # Extract the prey and predator solutions
    prey = y_vals[:, 0]
    predator = y_vals[:, 1]

    # Plot the results
    plt.figure(figsize=(8, 5))
    plt.plot(t_vals, prey, label="Prey", color="blue")
    plt.plot(t_vals, predator, label="Predator", color="red")
    plt.title("Predator-Prey Model (Explicit Euler)")
    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.legend()
    plt.grid(True)
    plt.show()
