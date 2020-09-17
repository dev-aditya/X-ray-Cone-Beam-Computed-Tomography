"""Required Imports"""
from numpy import array, zeros, transpose, squeeze

def shepp_logan():
	"""
	The lengths are in "cm" and angles i degrees
	"""
	xhat = array([0, 0, 0.22, - 0.22, 0, 0, 0, - 0.08, 0, 0.06]) * 10
	xhat = transpose(xhat)  ## you can also use numpy reshape function or method
	yhat = array([0, - 0.0184, 0, 0, 0.35, 0.1, - 0.1, - 0.605, - 0.605, - 0.605]) * 10
	yhat = transpose(yhat)
	zhat = array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) * 10
	zhat = transpose(zhat)
	a = array([0.69, 0.6624, 0.11, 0.16, 0.21, 0.046, 0.046, 0.046, 0.023, 0.023]) * 10
	a = transpose(a)
	b = array([0.92, 0.874, 0.31, 0.41, 0.25, 0.046, 0.046, 0.023, 0.023, 0.046]) * 10
	b = transpose(b)
	c = array([.9, .88, .21, .22, .35, .046, .02, .02, .1, .1]) * 10
	c = transpose(c)
	theta = squeeze(zeros((10, 1)))  ## squeeze -> to remove the dimension of length 1 form the array
	theta = transpose(theta)
	phi = array([0, 0, -18, 18, 0, 0, 0, 0, 0, 0])
	phi = transpose(phi)
	psi = squeeze(zeros((10, 1)))
	psi = transpose(psi)
	mu = array([1, - 0.8, - 0.2, - 0.2, - 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
	mu = transpose(mu)

	return array([xhat, yhat, zhat, a, b, c, theta, phi, psi, mu])


def yu_ye_wang():
	xhat = transpose(array([0, 0, -.22, .22, 0, 0, -.08, .06, 0.06, 0]))
	yhat = transpose(array([0, 0, 0, 0, 0.35, 0.1, -0.65, -0.65, -0.105, 0.1]))
	zhat = transpose(array([0, 0, -0.25, -0.25, -0.25, -0.25, -0.25, -0.25, 0.625, 0.625]))

	theta = transpose(squeeze(zeros((10, 1))))
	psi = transpose(squeeze(zeros((10, 1))))
	phi = transpose(array([0, 0, 108, 72, 0, 0, 0, 90, 90, 0]))

	a = transpose(array([0.69, 0.6624, 0.41, 0.31, 0.21, 0.046, 0.23, 0.23, 0.04, 0.56]))
	b = transpose(array([0.92, 0.874, 0.16, 0.11, 0.25, 0.046, 0.023, 0.023, 0.04, 0.56]))
	c = transpose(array([0.9, 0.88, 0.21, 0.22, 0.5, 0.046, 0.02, 0.02, 0.10, 0.10]))

	mu = transpose(array([1, -0.8, -.2, -.2, .2, .2, .1, .1, .2, -.2]))

	return array([xhat, yhat, zhat, a, b, c, theta, phi, psi, mu])



if __name__ == "__main__":
	import sys
	print('Python %s on %s' % (sys.version, sys.platform))
	print(f"Path to the file {__file__}")
# print(f"The shepp-logan phantom attributes: \n {shepp_logan()}")
