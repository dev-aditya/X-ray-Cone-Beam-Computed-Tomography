from numpy import pi, cos, sin, zeros, meshgrid, arange, array


########## PHANTON CONSTRUCTOR ############

def phantom_const(x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz, phantom_features):
	"""
	Parameters:
		x_cent : The center x coordinate of phantom in Nx*Ny*Nz volume geometry
		y_cent : The center y coordinate of phantom in Nx*Ny*Nz volume geometry
		z_cent : The center z coordinate of phantom in Nx*Ny*Nz volume geometry
		dx : voxel size in x direction
		dy : voxel size in y direction
		dz : voxel size in z direction
		phantom_att : The attributes of the phantom
		-----------------------------------------------------------------------
		Returns:
			Phantom of given shape, features and voxel size
	"""

	p = zeros((Nx, Ny, Nz))
	x = -x_cent * dx + arange(0, Nx) * dx
	y = -y_cent * dy + arange(0, Ny) * dy
	z = -z_cent * dz + arange(0, Nz) * dz

	x_coord, y_coord, z_coord = meshgrid(x, y, z)

	for k in range(len(phantom_features[0])):
		x0 = phantom_features[0][k]
		y0 = phantom_features[1][k]
		z0 = phantom_features[2][k]
		a = phantom_features[3][k]
		b = phantom_features[4][k]
		c = phantom_features[5][k]
		theta = phantom_features[6][k] * pi / 180
		phi = phantom_features[7][k] * pi / 180
		psi = phantom_features[8][k] * pi / 180
		f0 = phantom_features[9][k]

		e1 = array([cos(phi) * cos(theta), sin(phi) * cos(theta), -sin(theta)])
		e2 = array([-sin(phi), cos(phi), 0])

		n0 = cos(psi) * e1 + sin(psi) * e2
		n1 = -sin(psi) * e1 + cos(psi) * e2
		n2 = array([cos(phi) * sin(theta), sin(theta) * sin(phi), cos(theta)])

		aa = 1 / a ** 2
		bb = 1 / b ** 2
		cc = 1 / c ** 2

		p1 = n0[0] ** 2 * aa + n1[0] ** 2 * bb + n2[0] ** 2 * cc
		p2 = n0[1] ** 2 * aa + n1[1] ** 2 * bb + n2[1] ** 2 * cc
		p3 = n0[2] ** 2 * aa + n1[2] ** 2 * bb + n2[2] ** 2 * cc
		p4 = n0[0] * n0[1] * aa + n1[0] * n1[1] * bb + n2[0] * n2[1] * cc
		p5 = n0[0] * n0[2] * aa + n1[0] * n1[2] * bb + n2[0] * n2[2] * cc
		p6 = n0[2] * n0[1] * aa + n1[2] * n1[1] * bb + n2[2] * n2[1] * cc

		equation1 = p1 * (x_coord - x0) ** 2 + p2 * (y_coord - y0) ** 2 + p3 * (z_coord - z0) ** 2
		equation2 = p4 * (x_coord - x0) * (y_coord - y0) + p5 * (x_coord - x0) * (z_coord - z0) + p6 * (y_coord - y0) \
		            * (z_coord - z0)
		equation = equation1 + 2 * equation2

		ell = zeros(equation.shape)
		ell[equation < 1.0] = 1

		p = p + f0 * ell

	return p


if __name__ == "__main__":
	import phantom_features as features
	from numpy import flipud, squeeze
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	#%%%%%%%%%%%%%%%
	mpl.rcParams['figure.dpi'] = 600
	#%%%%%%%%%%%%%%%
	Nx = 256
	Ny = 256
	Nz = 128
	#%%%%%%%%%%%%%%%
	x_cent = int(Nx / 2)
	y_cent = int(Ny / 2)
	z_cent = int(Nz / 2)
	#%%%%%%%%%%%%%%%
	dx = 20 / Nx
	dy = 20 / Ny
	dz = 20 / Nz
	# %%%%%%%%%%%%%%%
	phantom_features = features.shepp_logan()
	phantom = phantom_const(x_cent, y_cent, z_cent, dx, dy, dz, Nx, Ny, Nz, \
	                        phantom_features)
	# %%%%%%%%%%%%%%%
	## squeeze to remove the dimension of length 1
	plt.imshow(flipud(squeeze(phantom[:, y_cent, :])), vmin = phantom.min(), vmax = phantom.max() \
	           , aspect = "equal")
	plt.show()
	plt.imshow(flipud(squeeze(phantom[x_cent, :, :])), vmin = phantom.min(), vmax = phantom.max() \
	           , aspect = "equal")
	plt.show()
	plt.imshow(flipud(squeeze(phantom[:, :, z_cent])), vmin = phantom.min(), vmax = phantom.max() \
	           , aspect = "equal")
	plt.show()
