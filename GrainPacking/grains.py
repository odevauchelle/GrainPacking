from pylab import *
from shapely.geometry import Polygon
from shapely.affinity import translate, rotate
from random import randrange

#############################
#
# Classes
#
#############################

class solid() :

    def __init__( self, **kwargs ) :
        self.Polygon = Polygon( **kwargs )

    def translate( self, **kwargs ) :
        self.Polygon = translate( self.Polygon, **kwargs )

    def rotate( self, **kwargs ) :
        self.Polygon = rotate( self.Polygon, **kwargs )

    def get_center( self ) :
        # return mean( self.Polygon.exterior.xy, axis = 1 )
        return self.Polygon.centroid.coords[0]

    def plot( self, ax = None, **kwargs ) :

        if ax is None :
            ax = gca()

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

def create_grain( npts = 8, radius = 1., center = (0.,0.), roughness = .2 ) :

    theta = sort( rand(npts)*2*pi )
    r = ( 1 + roughness*rand( npts ) )/( 1 + roughness/2 )
    shell = array(center).T + radius*array( [ cos(theta), sin(theta) ] ).T

    return solid( shell = shell )

def create_box( size = (1,1), thickness = .5, center = ( 0., 0. ) ) :

    corners = [ [ -1, -1 ],  [ 1, -1 ], [ 1, 1 ], [-1, 1] ]
    hole = [ array( point ).T*.5*array(size).T + array( center ).T for point in corners ]
    shell = [ array( point ).T*.5*( array(size).T*( 1 + thickness ) ) + array( center ).T for point in corners ]

    return solid( shell = shell, holes = [hole] )

def intersection_area( *solids ) :

    intersection = solids[0].Polygon

    for solid in solids[1:] :
        intersection = intersection.intersection( solid.Polygon )

    return intersection.area

def intersection_matrix( solids ) :

    n = len(solids)
    M = eye( n, n )*0.

    for i in range(n) :
        for j in range(i) :
            M[ i, j ] = intersection_area( solids[i], solids[j] )

    return M

def contact_Hamiltonian( solids ) :
    return sum( intersection_matrix( solids ) )

def gravity_Hamiltonian( solids, gravity = ( 0., -1. ) ) :
    return -sum( [ dot( the_solid.get_center(), gravity ) for the_solid in solids ] )

def Glauber_step( box, grains, Hamiltonian, beta, dx, dtheta = 0, current_energy = None ) :

    if current_energy is None :
        current_energy = Hamiltonian( box, grains )

    # pick a grain
    i = randrange( len(grains) )

    # choose a step
    dz = dx*rand()*exp( 1j*2*pi*rand() )

    # move the grain
    grains[i].translate( xoff = real(dz), yoff = imag(dz) )

    # rotate the grain
    if dtheta > 0 :
        grains[i].rotate( angle = dtheta*( rand() - .5 ), use_radians = True )

    # calculate the energy cost
    old_energy = current_energy
    current_energy = Hamiltonian( box, grains )
    dE = current_energy - old_energy

    # calculate probability
    p = exp( -beta*dE )/( 1 + exp( -beta*dE ) )

    if rand() < p : # keep
        return grains, current_energy

    else : # drop
        grains[i].translate( xoff = -real(dz), yoff = -imag(dz) )
        return grains, old_energy

def plot_grains( grains, ax = None, labels = False, **kwargs ) :

    if ax is None :
        ax = gca()

    for i, the_grain in enumerate( grains ) :
        the_grain.plot( ax = ax, **kwargs )

        if labels :
            text( *the_grain.get_center(), str(i), ha = 'center', va = 'center' )



############################
#
# Try it out
#
#############################

if __name__ == '__main__':

    ax_phys = gca()

    grain_radius = .1

    def Hamiltonian( box, grains ) :
        return gravity_Hamiltonian( grains ) + 1e1/grain_radius**2*contact_Hamiltonian( [box] + grains )

    box = create_box()

    grains = []

    for _ in range(15) :
        grains += [ create_grain( center = rand( 2 )-.5, radius = grain_radius ) ]

    current_energy = Hamiltonian( box, grains )

    plot_grains( grains, ax = ax_phys, color = 'tab:orange', alpha = .2 )

    E = []

    for beta in linspace( 10, 100, 100*len(grains) ) :
        dx = 2/beta
        dtheta = 0*dx/grain_radius
        E += [ current_energy ]
        grains, current_energy = Glauber_step( box, grains, Hamiltonian, beta = beta, dx = dx, dtheta = dtheta, current_energy = current_energy )

    box.plot(ax = ax_phys, color = 'grey')
    plot_grains( grains, ax = ax_phys, color = 'tab:brown' )

    ax_phys.axis('scaled')
    ax_phys.axis('off')

    figure()
    ax_e = gca()

    ax_e.plot( E )

    show()
