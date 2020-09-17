from numpy import sin, cos, zeros, arange, meshgrid, pi, ones, squeeze, array, flipud
from scipy.interpolate import RectBivariateSpline, interp2d

def back_project(g, u_cent, v_cent, du, dv, R0, D, x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz):
	Nv, Nu, Np = G.shape
	u_1 = (arange(0, Nu) - u_cent) * du
	v_1 = (arange(0, Nv) - v_cent) * dv
	[u, v] = meshgrid(u_1, v_1)

	dlambda = 2 * pi / Np
	lambda_ = dlambda * arange(0, Np)

	eu = array([-sin(lambda_), cos(lambda_), zeros(len(lambda_))])
	ev = array([zeros(len(lambda_)), zeros(len(lambda_)), ones(len(lambda_))])
	ew = array([cos(lambda_), sin(lambda_), zeros(len(lambda_))])

	xc = (arange(0, Nx) - x_cent) * dx
	yc = (arange(0, Ny) - y_cent) * dy
	zc = (arange(0, Nz) - z_cent) * dz
	[x, y] = meshgrid(xc, yc)

	back_proj = squeeze(zeros((Nx, Ny, Nz)))

	for i in arange(0, Nz, 1):
		z = zc[i]

		for j in arange(0, len(lambda_)):
			lam = lambda_[j]
			alamx = R0 * cos(lam) * ones(x.shape)
			alamy = R0 * sin(lam) * ones(y.shape)
			ustar_num = (x - alamx) * eu[j, 1] + (y - alamy) * eu[j, 2] + z * eu[j, 3]
			ustar_den = (x - alamx) * ew[j, 1] + (y - alamy) * ew[j, 2] + z * ew[j, 3]
			ustar = -D * (ustar_num / ustar_den)
			vstar_num = (x - alamx) * ev[j, 1] + (y - alamy) * ev[j, 2] + z * ev[j, 3]
			vstar_den = ustar_den
			vstar = -D * (vstar_num / vstar_den)
			gmod = interp2d(u, v, g[:, :, j])
			gmod = gmod(ustar,vstar)
			dd = ((alamx - x) * ew[j, 1] + (alamy - y) * ew[j, 2] + (0 - z) * ew[j, 3]) ** 2
			B0 = (R0 * D) / dd
			back_proj[:, :, i] = back_proj[:, :, j] + B0 * gmod * dlambda
			print("On i = %s and j = %s" % (i,j))
	back_proj = back_proj / 2
	return back_proj


if __name__ == '__main__':
	import phantom_features as features
	import projection as proj
	import matplotlib.pyplot as plt
	import scipy.ndimage as ndimg

	D = 200

	Nu = 256
	Nv = 256
	Np = 100  # number of projections
	u_cent = Nu / 2
	v_cent = Nv / 2
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

	al = 0.5

	G = proj.projection(u_cent, v_cent, du, dv, Nu, Nv, Np, R0, D, features.shepp_logan())

	B = back_project(G, u_cent, v_cent, du, dv, R0, D, x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz)
	del G
	plt.imshow(flipud(ndimg.laplace(squeeze(B[:, :, z_cent]))), vmin = 0, vmax = 0.001)
	plt.show()
	plt.imshow(flipud(ndimg.laplace(squeeze(B[:, y_cent, :]))), vmin = 0, vmax = 0.001)
	plt.show()
	plt.imshow(flipud(ndimg.laplace(squeeze(B[x_cent, :, :]))), vmin = 0, vmax = 0.001)
	plt.show()
