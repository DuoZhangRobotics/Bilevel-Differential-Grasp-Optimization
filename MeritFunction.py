import torch


class MeritFunction(object):
    def __init__(self, function, constraints_func, x0, dx0, pho=0.01, tol=1e-5, type='l1', last_eta=torch.tensor(1, dtype=torch.double)):
        """

        Parameters
        ----------
        function            objective function to be optimized
        eta                 the weights or coefficient of the l1 penalty term
        constraints_func    constraints h(x) = 0 and g(x) <= 0
        type                type of merit function including l1 penalty merit function and augmented lagrangian merit
                            function
        """
        self.function = function
        self.constraints_func = constraints_func
        self.x0 = x0
        self.dx0 = dx0
        self.pho = pho
        self.type = type
        self.tol = tol
        self.dfdx = torch.autograd.grad(self.function(x0), x0)[0] @ self.dx0
        constraints = self.constraints_func(x0)
        inequalities = torch.where(constraints <= 0, torch.tensor(0, dtype=torch.double), constraints)
        self.penalty_norm = torch.norm(inequalities, p=1)
        self.last_eta = last_eta
        self.eta = self._initialize_eta()
        self.directional_derivative = self.dfdx - self.eta * self.penalty_norm
        self.converged = (torch.abs(self.directional_derivative.detach()).numpy() <= self.tol)

    def merit_function(self, x):
        constraints = self.constraints_func(x)
        inequalities = torch.where(constraints <= 0, torch.tensor(0, dtype=torch.double), constraints)
        if self.type == 'l1':
            return self.function(x) + self.eta * torch.norm(inequalities, p=1)

    def _initialize_eta(self):
        if self.penalty_norm == 0:
            eta = torch.tensor(0, dtype=torch.double)
            return eta
        eta = self.dfdx / ((1 - self.pho) * self.penalty_norm)
        if eta < self.last_eta:
            eta = self.last_eta
        self.last_eta = eta
        return eta
