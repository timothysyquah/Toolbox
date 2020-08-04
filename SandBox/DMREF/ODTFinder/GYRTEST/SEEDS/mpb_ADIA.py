import math
import meep as mp
from meep import mpb

# Dielectric spheres in a diamond (fcc) lattice.  This file is used in
# the "Data Analysis Tutorial" section of the MPB manual.

sqrt_half = math.sqrt(0.5)
geometry_lattice = mp.Lattice(
    basis_size=mp.Vector3(sqrt_half, sqrt_half, sqrt_half),
    basis1=mp.Vector3(0, 1, 1),
    basis2=mp.Vector3(1, 0, 1),
    basis3=mp.Vector3(1, 1, 0)
)

# Corners of the irreducible Brillouin zone for the fcc lattice,
# in a canonical order:
vlist = [
    mp.Vector3(0, 0.5, 0.5),        # X
    mp.Vector3(0, 0.625, 0.375),    # U
    mp.Vector3(0, 0.5, 0),          # L
    mp.Vector3(0, 0, 0),            # Gamma
    mp.Vector3(0, 0.5, 0.5),        # X
    mp.Vector3(0.25, 0.75, 0.5),    # W
    mp.Vector3(0.375, 0.75, 0.375)  # K
]
tick_labs = ['X', 'U', 'L', '$\Gamma$', 'X','W','K']

k_points = mp.interpolate(4, vlist)

# define a couple of parameters (which we can set from the command_line)
#eps = 11.56  # the dielectric constant of the spheres
#r = 0.25  # the radius of the spheres

#diel = mp.Medium(epsilon=eps)

# A diamond lattice has two "atoms" per unit cell:
#geometry = [mp.Sphere(r, center=mp.Vector3(0.125, 0.125, 0.125), material=diel),
#            mp.Sphere(r, center=mp.Vector3(-0.125, -0.125, -0.125), material=diel)]

# (A simple fcc lattice would have only one sphere/object at the origin.)

resolution = 16  # use a 16x16x16 grid
mesh_size = 5
num_bands = 10

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    k_points=k_points,
#    geometry=geometry,
    resolution=resolution,
    num_bands=num_bands,
    mesh_size=mesh_size,
    epsilon_input_file = 'epsilon.h5'
)


def main():
    # run calculation, outputting electric_field energy density at the U point:
    #ms.run(mpb.output_at_kpoint(mp.Vector3(0, 0.625, 0.375), mpb.output_dpwr))
    ms.run()
    ## get epsilon
    #md = mpb.MPBData(rectify=True, periods=3, resolution=32)
    #eps = ms.get_epsilon()
    #converted_eps = md.convert(eps)
   

    import matplotlib
    matplotlib.use('Agg') # this allows plotting on braid
    import matplotlib.pyplot as plt
    
    # plot epsilon
    #fig, ax = plt.subplots()
    #im = ax.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
    #plt.axis('off')
    #plt.colorbar(im,label='$\epsilon (r)$')
    #plt.savefig('epsilon.png')

    fig, ax = plt.subplots()
    #kz_points = [k[0] for k in k_points]
    freqs = ms.all_freqs
    gaps = ms.gap_list
    x = range(len(freqs))
    # Plot bands
    # Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
    for xz, tmz, tez in zip(x, freqs, freqs):
        ax.scatter([xz]*len(tmz), tmz, color='blue')
    ax.plot(freqs, color='blue')
    ax.set_xlim([x[0], x[-1]])
    #ax.set_ylim([0,0.3])

    # Plot gaps
    for gap in gaps:
        if gap[0] > 1:
            ax.fill_between(x, gap[1], gap[2], color='blue', alpha=0.2)


    # Plot labels
    #ax.text(12, 0.04, 'TM bands', color='blue', size=15)
    #ax.text(13.05, 0.235, 'TE bands', color='red', size=15)

    points_in_between = (len(freqs) - len(vlist)) / (len(vlist)-1)
    tick_locs = [i*points_in_between+i for i in range(len(vlist))]
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs, size=16)

    ax.set_xlabel('Wave vector ' + '$ka/2\pi$', size=16)
    ax.set_ylabel('Frequency ' + '$\omega a / 2\pi c$', size=16)
    ax.grid(True)
    plt.tight_layout()
    plt.savefig('bands.png')


if __name__ == '__main__':
    main()
