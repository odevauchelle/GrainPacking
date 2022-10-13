from pylab import *
from matplotlib import gridspec
import json

sys.path.append('../')
import GrainPacking as GP

###################################
#
# parameters
#
###################################

p = dict(
    nb_grains = 10,
    grain = dict( npts = 7  )
    )

p['grain']['radius'] = .1
p['epsilon'] = 0.1*p['grain']['radius']**2
p['dx'] = 0.1*p['grain']['radius']
p['dtheta'] = p['dx']/p['grain']['radius']
p['beta_range'] = [ 10, 1000 ]
p['time_steps'] = 3000*p['nb_grains']
phi = pi/4
p['gravity'] = [ sin(phi), -cos(phi)  ]
p['grain']['roughness'] = .8

beta_list = linspace( *p['beta_range'], p['time_steps'] )

############################
#
# Create system
#
#############################


def Hamiltonian( box, grains ) :
    return GP.gravity_Hamiltonian( grains, gravity = p['gravity'] ) + 1/(p['epsilon'])*GP.contact_Hamiltonian( [box] + grains )

box = GP.create_box()

grains = []

for _ in range(p['nb_grains']) :
    grains += [ GP.create_grain( center = .5*( rand( 2 ) -.5 ), **p['grain'] ) ]

energy = Hamiltonian( box, grains )

###############################
#
# Prepare figure
#
###############################


fig = figure( figsize = (11,5), tight_layout = True )

gs = gridspec.GridSpec(2, 2)

ax_phys = fig.add_subplot( gs[ :, 0])
ax_e = fig.add_subplot( gs[ 1, 1])
ax_temp = fig.add_subplot( gs[ 0, 1], sharex = ax_e )

box.plot(ax = ax_phys, color = 'grey')
ax_phys.plot( *array([(0,0), p['gravity']]).T, '--k', alpha = .1 )
GP.plot_grains( grains, ax = ax_phys, color = 'tab:orange', alpha = .2 )

ax_phys.axis('scaled')
ax_phys.axis('off')
ax_phys.set_xticks([])
ax_phys.set_yticks([])
ax_e.set_ylabel('Energy')
ax_e.set_xlabel('Time')
ax_e.set_xlim( 0, p['time_steps']/len(grains) )

ax_temp.set_ylabel(r'Temperature $1/\beta$')
ax_temp.set_yscale('log')
ax_temp.plot( arange(p['time_steps'])/p['nb_grains'], 1/beta_list, 'tab:green'  )

E = []
e_graph, = ax_e.plot( arange(len(E))/len(grains), E )
g_graphs = GP.plot_grains( grains, ax = ax_phys, color = 'tab:brown' )

pause(.01)
show( block = False)

i_plot = 0

def refresh_plots( E, grains ) :

    e_graph.set_data( arange(len(E))/len(grains), E )
    ax_e.relim()
    ax_e.autoscale_view( scalex= False)

    for i, grain in enumerate(grains) :
        g_graphs[i].set_data( *grain.Polygon.exterior.xy )

    pause(0.01)

###############################
#
# Glauber loop
#
###############################


for beta in beta_list :

    i_plot += 1
    E += [ energy ]
    grains, energy = GP.Glauber_step( box, grains, Hamiltonian, beta = beta, dx = p['dx'], dtheta = p['dtheta'], energy = energy )

    if i_plot%300 == 0 :
        i_plot = 0
        refresh_plots( E, grains )


refresh_plots( E, grains )

###############################
#
# Save data
#
###############################

data_file = './grains.json'

data = dict( parameters = p )
data['box'] = box.dumps()
data['grains'] = [ grain.dumps() for grain in grains ]

if input('Save as ' + data_file + '?' ) in ( '', 'y' ) :

    with open( data_file, 'w' ) as the_file :
        json.dump( data, the_file )

    print('Saved.')

else :
    print('Not saved.')


# fig_phys.savefig('./triangles.svg', bbox_inches = 'tight')
