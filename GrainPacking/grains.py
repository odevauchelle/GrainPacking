import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mpl_Polygon
from random import randrange
from functools import reduce
from itertools import combinations
from shapely.geometry import Polygon as shp_Polygon
from shapely.affinity import translate, rotate
from shapely import wkt as shp_wkt
# from multiprocessing import Pool

#############################
#
# Classes
#
#############################

class solid() :

    def __init__( self, polygon = None, wkt_polygon = None, **kwargs ) :

        if polygon is None and wkt_polygon is None :
            self.Polygon = shp_Polygon( **kwargs )

        elif wkt_polygon is None :
            self.Polygon = polygon

        else : # use wkt_polygon
            self.Polygon = shp_wkt.loads( wkt_polygon )

        self.centroid = np.array( self.Polygon.centroid.coords[0] )

    def translate( self, **kwargs ) :
        self.Polygon = translate( self.Polygon, **kwargs )
        self.centroid[0] += kwargs['xoff']
        self.centroid[1] += kwargs['yoff']

    def rotate( self, **kwargs ) :
        self.Polygon = rotate( self.Polygon, origin = 'centroid', **kwargs )

    def get_center( self ) :
        # return mean( self.Polygon.exterior.xy, axis = 1 )
        return self.Polygon.centroid.coords[0]
        # return self.centroid

    def dumps( self ) :
        return shp_wkt.dumps( self.Polygon )
        # return self.Polygon.wkt.dumps() later shapely version ?

    def get_patch( self, **kwargs ) :
        return mpl_Polygon( np.array( self.Polygon.exterior.xy ).T )


    def plot( self, ax = None, **kwargs ) :

        if ax is None :
            ax = plt.gca()

        plots = [ ax.plot( *self.Polygon.exterior.xy, **kwargs ) ]

        kwargs.update( dict( color =  plots[0][0].get_color() ) )

        for hole in self.Polygon.interiors :
            plots += [ ax.plot( *hole.xy, **kwargs ) ]

        return plots

############################
#
# Functions
#
############################

def create_grain( npts = 3, radius = 1., center = (0.,0.), roughness = 0. ) :

    theta = np.sort( np.linspace( 0, 2*np.pi, npts + 1 )[:-1] + roughness*2*np.pi/npts*np.random.rand( npts ) )
    r = radius*( 1 + roughness*np.random.rand( npts ) )/( 1 + roughness/2 )
    shell = np.array(center).T + np.array( [ r*np.cos(theta), r*np.sin(theta) ] ).T

    the_grain = solid( shell = shell )
    the_grain.rotate( angle = np.random.rand()*360 )
    return the_grain

def create_box( size = ( 1., 1. ), thickness = .3, center = ( 0., 0. ) ) :

    corners = [ [ -1, -1 ],  [ 1, -1 ], [ 1, 1 ], [-1, 1] ]
    hole = [ np.array( point ).T*.5*np.array( size ).T + np.array( center ).T for point in corners ]
    shell = [ np.array( point ).T*.5*( np.array( size ).T*( 1 + 2*thickness ) ) + np.array( center ).T for point in corners ]

    return solid( shell = shell, holes = [ hole ] )

def intersection_area( solids ) :

    intersection = solids[0].Polygon

    for solid in solids[1:] :
        intersection = intersection.intersection( solid.Polygon )

    return intersection.area

def intersection_matrix( solids, M = None ) :

    if M is None :
        n = len(solids)
        M = eye( n, n )*0.

    for i in range(n) :
        for j in range(i) :
            M[ i, j ] = intersection_area( [ solids[i], solids[j] ] )

    return M

def contact_Hamiltonian( solids, multiprocessing_pool = None ) :

    if multiprocessing_pool is None :
        area = sum( map( intersection_area, combinations( solids, 2 ) ) )
    else :
        area = sum( multiprocessing_pool.map( intersection_area, combinations( solids, 2 ) ) )

    return area
    # return sum( intersection_matrix( solids ) )

def gravity_Hamiltonian( solids, gravity = ( 0., -1. ) ) :
    return -np.dot( np.sum( [ the_solid.centroid for the_solid in solids ], axis = 0 ), gravity )

def Glauber_step( box, grains, Hamiltonian, beta, dx, dtheta = 0., energy = None ) :

    if energy is None :
        energy = Hamiltonian( box, grains )

    dangle = None
    dz = None

    # pick a grain
    i = randrange( len(grains) )

    # rotate the grain...
    if dtheta > 0 :
        if np.random.choice( [True, False] ) :
            dangle = dtheta*( np.random.rand() - .5 )*180/np.pi
            grains[i].rotate( angle = dangle )

    # ...or move the grain
    if dangle is None :
        dz = dx*np.random.rand()*np.exp( 1j*2*np.pi*np.random.rand() )
        grains[i].translate( xoff = np.real(dz), yoff = np.imag(dz) )

    # calculate the energy np.cost
    old_energy = energy
    energy = Hamiltonian( box, grains )
    dE = energy - old_energy

    # calculate probability
    p = np.exp( -beta*dE )/( 1 + np.exp( -beta*dE ) )

    if np.random.rand() < p : # keep
        return grains, energy

    else : # drop

        try :
            grains[i].translate( xoff = -np.real(dz), yoff = -np.imag(dz) )

        except :
            grains[i].rotate( angle = -dangle )

        return grains, old_energy

def plot_grains( grains, ax = None, labels = False, **kwargs ) :

    if ax is None :
        ax = plt.gca()

    graphs = []

    for i, the_grain in enumerate( grains ) :
        graphs += [ the_grain.plot( ax = ax, **kwargs )[0][0] ]

        if labels :
            ax.text( *the_grain.centroid, str(i), ha = 'center', va = 'center' )

    return graphs



############################
#
# Try it out
#
#############################

if __name__ == '__main__':


    fig_phys = plt.figure()
    ax_phys = plt.gca()

    npts = 3
    grain_radius = 1/( 2*5*np.sin( 2*np.pi/npts )  )
    epsilon = .1

    def Hamiltonian( box, grains ) :
        return gravity_Hamiltonian( grains ) + 1/(epsilon*grain_radius**2)*contact_Hamiltonian( [box] + grains )

    box = create_box()

    grains = []


    for _ in range(10) :
        grains += [ create_grain( npts = npts, center = .5*( np.random.rand( 2 ) -.5 ), radius = grain_radius ) ]

    energy = Hamiltonian( box, grains )

    box.plot(ax = ax_phys, color = 'grey')
    plot_grains( grains, ax = ax_phys, color = 'tab:orange', alpha = .2 )

    ax_phys.axis('scaled')
    ax_phys.axis('off')
    ax_phys.set_xticks([])
    ax_phys.set_yticks([])

    plt.pause(.01)
    plt.show( block = False)
    E = []

    for beta in np.linspace( 30, 1000, 300*len(grains) ) :

        dx = 0.1*grain_radius
        dtheta = dx/grain_radius
        E += [ energy ]
        grains, energy = Glauber_step( box, grains, Hamiltonian, beta = beta, dx = dx, dtheta = dtheta, energy = energy )

    plot_grains( grains, ax = ax_phys, color = 'tab:brown' )


    plt.figure()

    ax_e = plt.gca()
    ax_e.plot( np.arange(len(E))/len(grains), E )
    ax_e.set_ylabel('Energy')
    ax_e.set_xlabel('time')

    # fig_phys.savefig('./triangles.svg', bbox_inches = 'tight')

    plt.show()
