from __future__ import division

import meep as mp
from meep import mpb

# Compute band structure for a square lattice of dielectric rods
# in air.

# Define various parameters with define_param so that they are
# settable from the command_line (with mpb <param>=<value>):
#r = 0.2  # radius of the rods
#eps = 11.56  # dielectric constant
k_interp = 4  # number of k points to interpolate

#GaAs = mp.Medium(epsilon=eps)

geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1))  # 2d cell

#geometry = [mp.Cylinder(r, material=GaAs)]

#Gamma = mp.Vector3()
#X = mp.Vector3(0.5, 0)
#M = mp.Vector3(0.5, 0.5)
vlist = [
    mp.Vector3(0.0, 0.0),        # Gamma
    mp.Vector3(0.5, 0.0),        # X
    mp.Vector3(0.5, 0.5),        # M
    mp.Vector3(0.0, 0.0)        # Gamma
]
tick_labs = ['$\Gamma$', 'X', 'M','$\Gamma$']

k_points = mp.interpolate(k_interp, vlist)

resolution = 32
num_bands = 8

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
#    geometry=geometry,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands,
    epsilon_input_file = 'epsilon.h5'
)


def main():
    ms.run_te()
    te_freqs = ms.all_freqs
    te_gaps = ms.gap_list
    ms.run_tm()
    tm_freqs = ms.all_freqs
    tm_gaps = ms.gap_list
    k_points = ms.k_points

    ms.display_eigensolver_stats()

    # start plot
    md = mpb.MPBData(rectify=True, periods=3, resolution=32)
    eps = ms.get_epsilon()
    converted_eps = md.convert(eps)

    import matplotlib
    matplotlib.use('Agg') # this allows plotting on braid
    import matplotlib.pyplot as plt
    # plot epsilon
    fig, ax = plt.subplots()
    im = ax.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
    plt.axis('off')
    plt.colorbar(im,label='$\epsilon (r)$')
    plt.savefig('epsilon.png')


    fig, ax = plt.subplots(figsize=(4,4))
    x = range(len(tm_freqs))
    #kz_points = [k[0] for k in k_points]
    # Plot bands
    # Scatter plot for multiple y values, see https://stackoverflow.com/a/34280815/2261298
    for xz, tmz, tez in zip(x, tm_freqs, te_freqs):
        ax.scatter([xz]*len(tmz), tmz, color='blue')
        #ax.scatter([xz]*len(tez), tez, color='red', facecolors='none')
    ax.plot(tm_freqs, color='blue')
    #ax.plot(te_freqs, color='red')
    #ax.set_ylim([0, 1])
    ax.set_xlim([x[0], x[-1]])
    #ax.set_ylim([0,0.3])

    # Plot gaps
    for gap in tm_gaps:
        if gap[0] > 1:
            ax.fill_between(x, gap[1], gap[2], color='blue', alpha=0.2)
    #for gap in te_gaps:
    #    if gap[0] > 1:
    #        ax.fill_between(x, gap[1], gap[2], color='red', alpha=0.2)

    # Plot labels
    #ax.text(12, 0.04, 'TM bands', color='blue', size=15)
    #ax.text(13.05, 0.235, 'TE bands', color='red', size=15)

    points_in_between = (len(tm_freqs) - len(vlist)) / (len(vlist)-1)
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
