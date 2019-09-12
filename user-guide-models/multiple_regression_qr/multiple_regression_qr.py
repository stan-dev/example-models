import pystan
import numpy as np

N = 100
K = 4
alpha = -3
beta = np.array([1,2,3,4], dtype="float")
sigma = 5

x = np.random.uniform(size = N * K)
x = x.reshape((N, K))
y = np.random.normal(size = N, loc=alpha + x.dot(beta), scale = sigma)

stan_data = {'N': N, 'K': K, 'x': x, 'y': y}

model = pystan.StanModel(file='multiple_regression_qr.stan')

fit = model.sampling(data=stan_data)
print(fit)
