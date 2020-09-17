from numpy import zeros, sin, cos, pi, abs, sqrt, sum, squeeze, arange


def ramp(proj, u_cent, v_cent, du, dv, src_obj, src_det):
	[nv, nu, nw] = proj.shape
	filtered = squeeze(zeros(proj.shape))
	g = squeeze(zeros((1, 256)))

	for deg in arange(0, nw):
		for i in arange(0, nu):
			for j in arange(0, nv):

				u = (i - 1 - u_cent) * du
				v = (j - 1 - v_cent) * dv

				for k in arange(0, nu):
					u_1 = du * (k - 1 - u_cent)
					dist = src_det / sqrt(u_1 ** 2 + v ** 2 + src_det ** 2) * du

					t_1 = pi * (u - u_1) / du
					t_2 = t_1 + pi
					t_3 = t_1 - pi
					if abs(t_1) < 1e-4:
						r_1 = 0.5
					elif abs(t_1) >= 10 ** (-4):
						r_1 = sin(t_1) / t_1 + (cos(t_1) - 1) / (t_1 ** 2)

					if abs(t_2) < 1e-4:
						r_2 = 0.5
					elif abs(t_2) >= 1e-4:
						r_2 = sin(t_2) / t_2 + (cos(t_2) - 1) / (t_2 ** 2)

					if abs(t_3) < 1e-4:
						r_3 = 0.5

					elif abs(t_3) >= 1e-4:
						r_3 = sin(t_3) / t_3 + (cos(t_3) - 1) / (t_3 ** 2)

					ram = src_obj * r_1 + 0.5 * (1 - src_obj) * (r_2 + r_3)
					ramp = 0.5 / du ** 2 * ram
					g[k] = ramp * dist * proj[j, k, deg]
			filtered[j, i, deg] = sum(g)

	return filtered
