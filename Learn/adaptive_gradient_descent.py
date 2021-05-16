# %%
import numpy as np
import scipy.optimize


def relative_grad(func, x, args, scale=0.001):
    """Calculate two-point approximation of gradient, with a relative step size."""
    step_mat = np.diag(x * scale)
    grad = np.zeros(len(x))
    for idx, steps in enumerate(step_mat):
        grad[idx] = (func(x + steps, *args) - func(x, *args)) / (steps[idx])

    return grad


# Armijo-Goldstein and positive
def condition_agp(func, grad, x, step_size, args, c=0.01, verbose=False):
    """Armijo-Goldstein condition for step size. c can be tweaked, but it's hard.
    Also checks for positive x.
    """

    # x must be positive after step
    if not all(x - grad * step_size > 0):
        if verbose:
            print(f"Negative x for: {step_size}")
        return False

    f_1 = func(x - step_size * grad, *args)
    f_0 = func(x, *args)

    # See https://en.wikipedia.org/wiki/Backtracking_line_search
    ag = f_1 <= f_0 - c * step_size * np.sum(grad ** 2)

    if not ag and verbose:
        print(f"Armijo-Goldstein failed for step size: {step_size}")

    return ag



def adaptive_agp_gradient_descent(func, x0, args, max_steps=100, step_size=1e-4, verbose=False):
    """Gradient descent with adaptive step size and positive constraint."""
    x = x0
    for idx in range(max_steps):
        grad = relative_grad(func, x, args, scale=0.0001)
        if condition_agp(func, grad, x, step_size, args, verbose=verbose):
            # Try increasing step size
            while condition_agp(func, grad, x, 2 * step_size, args, verbose=verbose):
                step_size = 2 * step_size

            x = x - step_size * grad
            if verbose:
                print(f"Function value: {func(x, *args)}")
                print(f"x: {x}, grad: {grad}")
        else:
            # Step size was to big! Half it until it passes
            step_size = step_size / 2
            while not condition_agp(func, grad, x, step_size, args, verbose=verbose):
                step_size = step_size / 2
            if verbose:
                print(f"Decreasing step to: {step_size}")

    return x, grad