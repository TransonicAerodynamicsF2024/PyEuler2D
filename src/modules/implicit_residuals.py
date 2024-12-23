import numpy as np

class ImplicitResidualSmoothing:
    def __init__(self, smoothing_coefficient=0.1, passes=1):
        """
        Initialize the IRS class.
        Args:
            smoothing_coefficient: Coefficient controlling smoothing intensity.
            passes: Number of smoothing passes.
        """
        self.epsilon = smoothing_coefficient
        self.passes = passes

    def thomas_algorithm(self, a, b, c, d):
        """
        Solve a tridiagonal system using the Thomas algorithm.
        """
        n = len(d)
        c_prime = np.zeros(n-1)
        d_prime = np.zeros(n)

        # Forward sweep
        c_prime[0] = c[0] / b[0]
        d_prime[0] = d[0] / b[0]
        for i in range(1, n-1):
            denominator = b[i] - a[i-1] * c_prime[i-1]
            c_prime[i] = c[i] / denominator
            d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / denominator
        d_prime[-1] = (d[-1] - a[-2] * d_prime[-2]) / (b[-1] - a[-2] * c_prime[-2])

        # Backward substitution
        x = np.zeros(n)
        x[-1] = d_prime[-1]
        for i in range(n-2, -1, -1):
            x[i] = d_prime[i] - c_prime[i] * x[i+1]

        return x

    def smooth(self, residual, mesh):
        """
        Apply IRS to the residuals in both i- and j-directions.
        """
        ni, nj = mesh.ni + 2 * mesh.n_ghosts, mesh.nj + 2 * mesh.n_ghosts
        smoothed_residual = np.copy(residual)

        for _ in range(self.passes):
            # Smoothing in i-direction
            for j in range(mesh.n_ghosts, nj - mesh.n_ghosts):
                a = -self.epsilon * np.ones(ni - 2 * mesh.n_ghosts - 1)
                b = (1 + 2 * self.epsilon) * np.ones(ni - 2 * mesh.n_ghosts)
                c = -self.epsilon * np.ones(ni - 2 * mesh.n_ghosts - 1)
                d = smoothed_residual[mesh.n_ghosts:ni-mesh.n_ghosts, j]
                smoothed_residual[mesh.n_ghosts:ni-mesh.n_ghosts, j] = self.thomas_algorithm(a, b, c, d)

            # Smoothing in j-direction
            for i in range(mesh.n_ghosts, ni - mesh.n_ghosts):
                a = -self.epsilon * np.ones(nj - 2 * mesh.n_ghosts - 1)
                b = (1 + 2 * self.epsilon) * np.ones(nj - 2 * mesh.n_ghosts)
                c = -self.epsilon * np.ones(nj - 2 * mesh.n_ghosts - 1)
                d = smoothed_residual[i, mesh.n_ghosts:nj-mesh.n_ghosts]
                smoothed_residual[i, mesh.n_ghosts:nj-mesh.n_ghosts] = self.thomas_algorithm(a, b, c, d)

        return smoothed_residual

