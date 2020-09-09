import torch


class MeritFunction(object):
    def __init__(self, function, eta, constraints_func, type='l1'):
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
        self.eta = eta
        self.constraints_func = constraints_func
        self.type = type

    def merit_function(self, x):
        if self.type == 'l1':
            return self.function(x) + self.eta * torch.norm(self.constraints_func(x), p=1)
        if self.type == 'alm':
            return self.function(x) + 0.5 * self.eta * torch.norm(self.constraints_func(x), p=2)

    def get_derivative(self, x):
        if self.type == 'l1':
            return torch.autograd.grad(self.function, x) + self.eta * torch.norm(self.constraints_func(x), p=1)
        if self.type == 'alm':
            return torch.autograd.grad(self.merit_function(x), x)

