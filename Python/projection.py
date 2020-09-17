from numpy import sin, cos, zeros, arange, meshgrid, pi, array, squeeze, matmul, transpose, real, flipud
from scipy import sqrt



def projection(u_cent, v_cent, du, dv, Nu, Nv, Np, R0, D, phant_feat):
	proj = squeeze(zeros((Nv, Nu, Np)))
	u_1 = (arange(0, Nu) - u_cent) * du
	v_1 = (arange(0, Nv) - v_cent) * dv

	[u, v] = meshgrid(u_1, v_1)

	dlambda = 2 * pi / Np
	lambda_ = dlambda * arange(0, Np,1)

	for i in arange(0, len(lambda_)):
		lam = lambda_[i]

		for k in arange(0, len(phant_feat[0])):
			x0 = phant_feat[0][k]
			y0 = phant_feat[1][k]
			z0 = phant_feat[2][k]
			a = phant_feat[3][k]
			b = phant_feat[4][k]
			c = phant_feat[5][k]
			theta = phant_feat[6][k] * pi / 180
			phi = phant_feat[7][k] * pi / 180
			psi = phant_feat[8][k] * pi / 180
			f0 = phant_feat[9][k]

			e1 = array([cos(phi) * cos(theta), sin(phi) * cos(theta), - sin(theta)])
			e2 = array([-sin(phi), cos(phi), 0])

			alpha = squeeze(zeros((Nv, Nu, 3)))
			n1 = cos(psi) * e1 + sin(psi) * e2
			n2 = -sin(psi) * e1 + cos(psi) * e2
			n3 = array([cos(phi) * sin(theta), sin(theta) * sin(phi), cos(theta)])

			a = 1 / a ** 2
			b = 1 / b ** 2
			c = 1 / c ** 2

			A = array([[a, 0, 0], [0, b, 0], [0, 0, c]])

			Q = array([n1, n2, n3])

			E0 = matmul(matmul(transpose(Q), A), Q)

			ew = array([cos(lam), sin(lam), 0])
			eu = array([-sin(lam), cos(lam), 0])
			ev = array([0, 0, 1])

			for j in range(0, 3):
				U = u * eu[j] + v * ev[j] - D * ew[j]
				t = u ** 2 + v ** 2 + D ** 2
				alpha[:, :, j] = t**(-1) * U
			alpha1 = alpha[:, :, 0]
			alpha2 = alpha[:, :, 1]
			alpha3 = alpha[:, :, 2]

			ahat = array([R0 * cos(lam), R0 * sin(lam), 0])
			xhat = ahat - array([x0, y0, z0])

			el0 = alpha1 * (alpha1 * E0[0, 0] + alpha2 * E0[0, 1] + alpha3 * E0[0, 2]) + \
			      alpha2 * (alpha1 * E0[1, 0] + alpha2 * E0[1, 1] + alpha3 * E0[1, 2]) + \
			      alpha3 * (alpha1 * E0[2, 0] + alpha2 * E0[2, 1] + alpha3 * E0[2, 2])

			el1 = alpha1 * (xhat[0] * E0[0, 0] + xhat[1] * E0[0, 1] + xhat[2] * E0[0, 2]) + \
			      alpha2 * (xhat[0] * E0[1, 0] + xhat[1] * E0[1, 1] + xhat[2] * E0[1, 2]) + \
			      alpha3 * (xhat[0] * E0[2, 0] + xhat[1] * E0[2, 1] + xhat[2] * E0[2, 2])

			el2 = matmul(matmul(xhat, E0), transpose(xhat)) - 1

			el = 1 / el0

			tp = (-el1 - real(sqrt(el1 ** 2 - el0 * el2))) * el
			tq = (-el1 + real(sqrt(el1 ** 2 - el0 * el2)))* el

			proj[:, :, k] = proj[:, :, k] + f0 * (tq - tp)

	return proj


if __name__ == "__main__":
	import phantom_features as features
	phant_feat = features.shepp_logan()
	D = 200

	Nu = 256
	Nv = 256
	Np = 100  # number of projections
	u_cent = int(Nu / 2)
	v_cent = int(Nv / 2)
	R0 = 200
	du = 40 / Nu
	dv = du

	Nx = 25
	Ny = 250
	Nz = 128
	x_cent = int(Nx / 2)
	y_cent = int(Ny / 2)
	z_cent = int(Nz / 2)

	dx = 20 / Nx
	dy = 20 / Nx
	dz = 20 / Nz

	import scipy.io as io
	proj = projection(u_cent, v_cent, du, dv, Nu, Nv, Np, R0, D, phant_feat)

	io.savemat("projection.mat", {"proj": proj})
	#import matplotlib.pyplot as plt

	#shape = proj.shape

	#plt.imshow(flipud(squeeze(proj[int(shape[0] / 2),:,:])),vmax = 0.2,vmin =0,aspect = "auto")
	#plt.show()
	#plt.imshow(squeeze(proj[:,int(shape[1] / 2),:]),vmax = 0.4,vmin =0,aspect ="auto")
	#plt.show()
	#plt.imshow(squeeze(proj[:,:,int(shape[2]/2)]),vmax = 0.01,vmin =0,aspect = "auto")
	#plt.show()