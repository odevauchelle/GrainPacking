import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mpl_Polygon
from random import randrange
# from functools import reduce
# from itertools import combinations
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

class system() :

    def __init__( self, box = None, grains = None, diff_Hamiltonian = None ) :

        self.box = box
        self.grains = grains
        self.set_dH( diff_Hamiltonian )

    def set_dH( self, diff_Hamiltonian ) :

        if diff_Hamiltonian is None :
            def dH(i) :
                return 0
        else :
            def dH( i ) :
                return diff_Hamiltonian( self.box, self.grains, i )

        self.dH = dH

    def set_generic_dH( self, epsilon = None, gravity = None ) :

        if epsilon is None and gravity is None :
            def dH(i) :
                return 0

        elif gravity is None :
            def dH( i ) :
                return 1/epsilon*diff_contact_Hamiltonian( [self.box] + self.grains, i + 1 )

        elif epsilon is None :
            def dH( i ) :
                return diff_gravity_Hamiltonian( self.grains, i, gravity )

        else :
            def dH( i ) :
                return diff_gravity_Hamiltonian( self.grains, i, gravity ) + 1/epsilon*diff_contact_Hamiltonian( [self.box] + self.grains, i + 1 )


        self.dH = dH

    def get_energy( self ) :
        return sum( map( self.dH, range( len( self.grains ) ) ) )

    def plot_grains( self, **kwargs ) :
        return plot_grains( self.grains, **kwargs )

    def plot_box( self, **kwargs ) :
        return self.box.plot( **kwargs )

    def Glauber_step( self, beta, dx, dtheta = 0. ) :

        dangle = None
        dz = None

        # pick a grain
        i = randrange( len(  self.grains ) )

        # calculate its local Energy
        dE = self.dH( i )

        # rotate the grain...
        if dtheta > 0 :
            if np.random.choice( [True, False] ) :
                dangle = dtheta*( np.random.rand() - .5 )*180/np.pi
                self.grains[i].rotate( angle = dangle )

        # ...or move the grain
        if dangle is None :
            dz = dx*np.random.rand()*np.exp( 1j*2*np.pi*np.random.rand() )
            self.grains[i].translate( xoff = np.real(dz), yoff = np.imag(dz) )

        # calculate the energy cost
        dE -= self.dH( i )
        dE *= -1

        # calculate probability
        p = np.exp( -beta*dE )/( 1 + np.exp( -beta*dE ) )

        if np.random.rand() > p : # drop, otherwise keep

            try :
                self.grains[i].translate( xoff = -np.real(dz), yoff = -np.imag(dz) )

            except :
                self.grains[i].rotate( angle = -dangle )

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

def intersection_matrix_line( solids, line_index ) :
    return [ intersection_area( [ solids[line_index], solid ] ) for solid in solids ]

def diff_contact_Hamiltonian( solids, index ) :
    return sum( intersection_matrix_line( solids, index ) )

def diff_gravity_Hamiltonian( solids, index, gravity = ( 0., -1. ) ) :
    return -np.dot( solids[index].centroid, gravity )

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
    epsilon = .1*grain_radius**2
    dx = 0.1*grain_radius
    dtheta = dx/grain_radius

    the_system = system()
    the_system.box = create_box()

    the_system.grains = []
    for _ in range(10) :
        the_system.grains += [ create_grain( npts = npts, center = .4*( np.random.rand( 2 ) -.5 ), radius = grain_radius ) ]

    the_system.set_generic_dH( epsilon = epsilon, gravity = [0,-1] )

    the_system.plot_box( ax = ax_phys, color = 'grey' )
    the_system.plot_grains( ax = ax_phys, color = 'tab:orange', alpha = .2 )

    ax_phys.axis('scaled')
    ax_phys.axis('off')
    ax_phys.set_xticks([])
    ax_phys.set_yticks([])

    plt.pause(.01)
    plt.show( block = False)

    for beta in np.linspace( 10, 1000, 300*len(the_system.grains) ) :

        the_system.Glauber_step(beta = beta, dx = dx, dtheta = dtheta )

    the_system.plot_grains( ax = ax_phys, color = 'tab:brown' )


    # fig_phys.savefig('./triangles.svg', bbox_inches = 'tight')

    plt.show()
