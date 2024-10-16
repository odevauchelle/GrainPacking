from pylab import *
from matplotlib import gridspec
import json
# from multiprocessing import Pool

sys.path.append('../../../GrainPacking/')
import GrainPacking as GP

####################################
# 0. Devauchelle & H. Chauvet
# October 2024
# Inspired by Sharon Glotzer's 2019 APS talk
# https://youtu.be/JW1L_vZ6K1M?feature=shared&t=1376
#
#####################################

###################################
#
# parameters
#
###################################

p = dict(
    nb_grains = 100,
    grain = dict( npts = 6 )
    )

p['grain']['radius'] = .07/sqrt(2)
p['epsilon'] = 0.1*p['grain']['radius']**2

p['dx'] = 2*p['grain']['radius'] # sets diffusivity
p['dtheta'] = 0.3*2*pi/p['grain']['npts']

p['beta_range'] = array([100,1000])
p['time_steps'] = 3000*p['nb_grains']
phi = 0
# p['gravity'] = [ sin(phi), -cos(phi)  ]
p['grain'].update( roughness = 0 )

beta_list = linspace( *p['beta_range'], p['time_steps'] )

############################
#
# Create system
#
#############################

s = GP.system()
s.box = GP.create_box()
s.grains = []

for _ in range(p['nb_grains']) :
    s.grains += [ GP.create_grain( center = .9*( rand( 2 ) -.5 ), **p['grain'] ) ]


def diff_Hamiltonian( box, grains, index ) : # only steric repulsion
    return 1/p['epsilon']*GP.diff_contact_Hamiltonian( [box] + grains, index + 1 )

s.set_dH( diff_Hamiltonian )

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

s.plot_box(ax = ax_phys, color = 'grey')
# ax_phys.plot( *array([(0,0), p['gravity']]).T, '--k', alpha = .1 )
# s.plot_grains( ax = ax_phys, color = 'tab:orange', alpha = .2 )

ax_phys.axis('scaled')
ax_phys.axis('off')
ax_phys.set_xticks([])
ax_phys.set_yticks([])
ax_e.set_ylabel('Energy')
ax_e.set_xlabel('Time')
ax_e.set_xlim( 0, p['time_steps']/len(s.grains) )

ax_temp.set_ylabel(r'Temperature $1/\beta$')
ax_temp.set_yscale('log')
ax_temp.plot( arange(p['time_steps'])/p['nb_grains'], 1/beta_list, 'tab:green'  )

E = []
time_E = []
e_graph, = ax_e.plot( time_E, E )
g_graphs = s.plot_grains( ax = ax_phys, color = 'tab:brown' )

pause(.01)
show( block = False)

def refresh_plots( time_E, E, grains ) :

    e_graph.set_data( time_E, E )
    ax_e.relim()
    ax_e.autoscale_view( scalex = False  )

    for i, grain in enumerate(grains) :
        g_graphs[i].set_data( *grain.Polygon.exterior.xy )

    pause(0.01)

###############################
#
# Glauber loop
#
###############################

i = 0

for beta in beta_list :

    i += 1
    s.Glauber_step( beta = beta, dx = p['dx'], dtheta = p['dtheta'] )

    if i%1000 == 0 :

        E += [ s.get_energy() ]
        time_E += [ i/len(s.grains) ]
        refresh_plots( time_E, E, s.grains )

E += [ s.get_energy() ]
time_E += [ i ]

refresh_plots( time_E, E, s.grains )

show(block = True)

###############################
#
# Save data
#
# ###############################

# data_file = './grains.json'

# data = dict( parameters = p )
# data['box'] = s.box.dumps()
# data['grains'] = [ grain.dumps() for grain in s.grains ]

# if input('Save as ' + data_file + '?' ) in ( '', 'y' ) :

#     with open( data_file, 'w' ) as the_file :
#         json.dump( data, the_file )

#     print('Saved.')

# else :
#     print('Not saved.')


# fig_phys.savefig('./triangles.svg', bbox_inches = 'tight')
