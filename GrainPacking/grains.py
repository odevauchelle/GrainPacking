from pylab import *
from shapely.geometry import Polygon
from shapely.affinity import translate, rotate
from shapely import wkt
from random import randrange

#############################
#
# Classes
#
#############################

class solid() :

    def __init__( self, polygon = None, wkt_polygon = None, **kwargs ) :

        if polygon is None and wkt_polygon is None :
            self.Polygon = Polygon( **kwargs )

        elif wkt_polygon is None :
            self.Polygon = polygon

        else : # use wkt_polygon
            self.Polygon = wkt.loads( wkt_polygon )

    def translate( self, **kwargs ) :
        self.Polygon = translate( self.Polygon, **kwargs )

    def rotate( self, **kwargs ) :
        self.Polygon = rotate( self.Polygon, origin = 'centroid', **kwargs )

    def get_center( self ) :
        # return mean( self.Polygon.exterior.xy, axis = 1 )
        return self.Polygon.centroid.coords[0]

    def dumps( self ) :
        return wkt.dumps( self.Polygon )
        # return self.Polygon.wkt.dumps() later shapely version ?

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
#
# def solid_loads( wkt_polygon ) :
#     return solid( polygon = wkt.loads( wkt_polygon ) )

def create_grain( npts = 3, radius = 1., center = (0.,0.), roughness = 0. ) :

    theta = sort( linspace( 0, 2*pi, npts + 1 )[:-1] + roughness*2*pi/npts*rand( npts ) )
    r = radius*( 1 + roughness*rand( npts ) )/( 1 + roughness/2 )
    shell = array(center).T + array( [ r*cos(theta), r*sin(theta) ] ).T

    the_grain = solid( shell = shell )
    the_grain.rotate( angle = rand()*360 )
    return the_grain

def create_box( size = (1,1), thickness = .3, center = ( 0., 0. ) ) :

    corners = [ [ -1, -1 ],  [ 1, -1 ], [ 1, 1 ], [-1, 1] ]
    hole = [ array( point ).T*.5*array(size).T + array( center ).T for point in corners ]
    shell = [ array( point ).T*.5*( array(size).T*( 1 + 2*thickness ) ) + array( center ).T for point in corners ]

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

def Glauber_step( box, grains, Hamiltonian, beta, dx, dtheta = 0, energy = None ) :

    if energy is None :
        energy = Hamiltonian( box, grains )

    dangle = None
    dz = None

    # pick a grain
    i = randrange( len(grains) )

    # rotate the grain...
    if dtheta > 0 :
        if rand() > .5 :
            dangle = dtheta*( rand() - .5 )*180/pi
            grains[i].rotate( angle = dangle )

    # ...or move the grain
    if dangle is None :
        dz = dx*rand()*exp( 1j*2*pi*rand() )
        grains[i].translate( xoff = real(dz), yoff = imag(dz) )

    # calculate the energy cost
    old_energy = energy
    energy = Hamiltonian( box, grains )
    dE = energy - old_energy

    # calculate probability
    p = exp( -beta*dE )/( 1 + exp( -beta*dE ) )

    if rand() < p : # keep
        return grains, energy

    else : # drop

        try :
            grains[i].translate( xoff = -real(dz), yoff = -imag(dz) )

        except :
            grains[i].rotate( angle = -dangle )

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

    fig_phys = figure()
    ax_phys = gca()

    npts = 3
    grain_radius = 1/( 2*4*sin( 2*pi/npts )  )
    epsilon = .1

    def Hamiltonian( box, grains ) :
        return gravity_Hamiltonian( grains ) + 1/(epsilon*grain_radius**2)*contact_Hamiltonian( [box] + grains )

    box = create_box()

    grains = []


    for _ in range(7) :
        grains += [ create_grain( npts = npts, center = .5*( rand( 2 ) -.5 ), radius = grain_radius ) ]

    energy = Hamiltonian( box, grains )

    box.plot(ax = ax_phys, color = 'grey')
    plot_grains( grains, ax = ax_phys, color = 'tab:orange', alpha = .2 )

    ax_phys.axis('scaled')
    ax_phys.axis('off')
    ax_phys.set_xticks([])
    ax_phys.set_yticks([])

    pause(.01)
    show( block = False)
    E = []

    for beta in linspace( 5, 100, 4000*len(grains) ) :
        dx = 0.1*grain_radius#2/beta
        dtheta = dx/grain_radius
        E += [ energy ]
        grains, energy = Glauber_step( box, grains, Hamiltonian, beta = beta, dx = dx, dtheta = dtheta, energy = energy )

    plot_grains( grains, ax = ax_phys, color = 'tab:brown' )


    figure()

    ax_e = gca()
    ax_e.plot( arange(len(E))/len(grains), E )
    ax_e.set_ylabel('Energy')
    ax_e.set_xlabel('time')

    fig_phys.savefig('./triangles.svg', bbox_inches = 'tight')

    show()
