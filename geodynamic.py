#!/usr/local/bin/python
# Time-stamp: <2016-02-23 09:44:14 marine>
# Project : From geodynamic to Seismic observations in the Earth's inner core
# Author : Marine Lasbleis



import numpy as np
import matplotlib.pyplot as plt #for figures
from mpl_toolkits.basemap import Basemap #to render maps
import math
from scipy.integrate import ode
from scipy.optimize import fsolve

#personal routines
import positions


RICB = 1221.


def exact_translation(point, velocity, direction=positions.CartesianPoint(1,0,0)):
    x_0, y_0, z_0 = point.x, point.y, point.z
    mean_direction = np.sqrt(direction.x**2+ direction.y**2+ direction.z**2)
    a, b, c = direction.x/mean_direction, direction.y/mean_direction, direction.z/mean_direction

    solution_1 = x_0*a+y_0*b+z_0*c + np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**2))
    solution_2 = x_0*a+y_0*b+z_0*c - np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**    2))
    ## TO DO : verify that we can remove solution_2 ? Is solution_1 always max?
    return max(solution_1, solution_2)

def velocity(t, r, vt, omega):
    """ velocity, in cartesian coordinates 
    
    r is a cartesian position [x, y, z]
    vt is the translation velocity (numpy array (1,3))
    rotation velocity is $\omega \vec{e}_z \times \vec{r}$.
    """
    return vt + np.array([-omega * r[1], omega * r[0], 0.] )

def integration_trajectory(r0, t0, t1, vt, omega):
    """ 
    
    r0: initial position
    t0: initial time
    t1: tmax of the integration
    """

    r = ode(velocity).set_integrator('zvode', method='bdf')
    r.set_initial_value(r0, t0).set_f_params(vt, omega) # .set_f_params() if the function has any parameters
    return np.real(r.integrate(r.t+(t1-t0)))

def trajectory_r(r0, t0, t1, vt, omega):
    trajectory = integration_trajectory(r0, t0, t1, vt, omega)
    print "cart", trajectory[0], trajectory[1], trajectory[2]
    r, t, p = positions.from_cartesian_to_seismo(trajectory[0], trajectory[1], trajectory[2])
    print "seismo", r,t,p
    return r

def findIntersection(fun1,fun2,x0):
     return fsolve(lambda x : fun1(x) - fun2(x),x0)

def radius_ic(t):
    return rICB

def find_age(r0, t0, vt, omega):
    return fsolve(lambda x : trajectory_r(r0, t0, x, vt, omega)-radius_ic(x), 0.)


    

def numberpoints_surface(r, density):
    """ return the number of points needed to cover the surface at radius 'r' with a density 'density'
    
    r: radius (radius = 1 correspond to IC now)
    density: density of points at the surface (corresponding to r=1)
    """
    return int(4.*np.pi*r**2.*density)


class Model():
    """ Geodynamical modelling of a flow. """

    def __init__(self):
        """ Point is a positions.Point """
        self.points = []
        self.proxy = []
        self.old = []
        self.new = [] 
        self.time = None
        self.history = []
        self.radius = None
        self.Npoints = None
        self.density = None
        self.tmax, self.dt = None, None

    def initialisation(self, t0, R0, density):
        self.time = t0
        self.radius = R0
        self.density = density
        self.fill_sphere(density) #fill the initial sphere with points. If the initial sphere is too small, only the point in r=0 is set up.
        self.Npoints = len(self.points)

    def run(self, tmax, dt):
        t0 = self.time

        print " ================="
        print " ===== RUN ======="
        print " ================="
        print "from t=", t0, " to t=", tmax, ", with dt=", dt, "."
        for t in np.arange(t0, tmax, dt):
            print "t = ", t, ". ", self.Npoints, " points."
            self.old = self.points
            self.time = t
            self.advect()
            self.grow()
            self.check_inside()
            #self.add_points(self.density_surface)
            self.Npoints = len(self.points)
            assert(self.Npoints == len(self.proxy))

    def advect(self):
        """ advection of the point. Need to be implemented in the derived classes """
        raise NotImplementedError("need to implement advection() in derived class!")

    def grow(self):
        """ grow the radius of the core """
        raise NotImplementedError("need to implement grow() in derived class!")

    def check_inside(self):
        """ check if the points are inside or outside the core"""
        # find points outside
        index_outside = [] #store the index of points to remove
        for i, point in enumerate(self.points):
            if point.r > self.radius:
                index_outside.append(i)




        # remove them (store?)
        self.points = np.delete(self.points, index_outside)
        self.proxy = np.delete(self.proxy, index_outside)

    def extract_rtp(self):
        r, theta, phi = np.empty([self.Npoints,1]), np.empty([self.Npoints,1]), np.empty([self.Npoints,1])
        for i, point in enumerate(self.points):
            r[i] = point.r
            theta[i] = point.theta
            phi[i] = point.phi
        return r, theta, phi

    def extract_xyz(self):
        x, y, z = np.empty([self.Npoints,1]), np.empty([self.Npoints,1]), np.empty([self.Npoints,1])
        for i, point in enumerate(self.points):
            x[i] = point.x
            y[i] = point.y
            z[i] = point.z
        return x, y, z 

    def add_points(self, density_surface, value=0.):
        """ add points where there are few of them and on the surface if ic is growing """
        # on the surface
        S = 4. * np.pi * self.radius**2 # size of the surface
        n = int(np.ceil(S*density_surface)) #number of points to add on the surface. 
        for i in range(n):
            u, v = np.random.uniform(0, 1), np.random.uniform(0, 1)
            r = self.radius
            theta = 90. - 180./np.pi*np.arccos(2*v-1) #here is the latitude
            phi = 360. * u
            self.points = np.append(self.points, positions.SeismoPoint(r, theta, phi))
            self.proxy = np.append(self.proxy, value)
        self.Npoints = len(self.points)


    def fill_sphere(self, density, initial_value=0.):
        volume = 4./3. * np.pi *self.radius**3. #volume of the sphere
        n = int(math.ceil(density*volume))
        if n == 1:
            self.points = np.append(self.points, positions.SeismoPoint(0.,0.,0.)) #central point
        else:
            # to get a uniform repartition, make it uniform in the (x, y, z) space, but with more points that needed, and then remove all the extra points. Check if at least one point is in the space, and if not add a point at the center. 
            n = int(math.ceil(density*volume*6./np.pi))
            for i in range(n):
                P = positions.CartesianPoint(self.radius*np.random.uniform(-1., 1.), self.radius*np.random.uniform(-1., 1.), self.radius*np.random.uniform(-1., 1.))
                if P.r <= self.radius:
                    self.points = np.append(self.points, P)
            self.Npoints = len(self.points)
            if self.Npoints <= 1:
                self.points = [positions.SeismoPoint(0.,0.,0.)]
                self.Npoints = 1
        self.proxy = initial_value* np.ones_like(self.points)                


    #def translation(self, velocity, dt):
    #    """ translation of the point.
    #
    #         velocity: tupple with 3 components (cartesian coordinates. x,y,z)
    #       """
    #       self.new.x = self.old.x + velocity*dt

class Translation(Model):
    ## Not up-to-date
    ### the __init__ method is directly obtained from the class Point_evolution, so no need to redefined it in a derived class
    ## def __init__(self, Point, time)
    ##     Point_evolution.__init__(self, Point, time)


    def initialisation(self, t0, R0, density):
        Model.initialisation(self, t0, R0, density)
        self.twoD_equatorial_plot()

    def advect(self):
        """ advection of the point. """

        dx = self.translation(self.dt, self.vt, self.direction)
        for i, point in enumerate(self.points): 
            #self.points[i].move(dx)
            self.points[i].move(dx)

    def grow(self):
        pass

    def translation(self, dt, velocity, direction):
        direction = direction / np.sqrt(direction[0]**2+direction[1]**2+direction[2]**2) #check if the direction has a size of 1
        dx = np.array([velocity*direction[0], velocity*direction[1], velocity*direction[2]])*dt
        return dx

    def run(self, tmax, dt):
        t0 = self.time
        fig, ax = plt.subplots()
        plt.axis('equal')
        print " ================="
        print " ===== RUN ======="
        print " ================="
        print "from t=", t0, " to t=", tmax, ", with dt=", dt, "."
        for t in np.arange(t0, tmax, dt):
            print "t = ", t, ". ", self.Npoints, " points."
            self.old = self.points
            self.time = t
            self.advect()
            self.check_inside()
            self.Npoints = len(self.points)
            assert(self.Npoints == len(self.proxy))
            x, y, z = self.extract_xyz()
            ax.plot(x, y, '.')

    def twoD_equatorial_plot(self, *args):
        if len(args) == 0: #initialisation of the plot
            fig, ax = plt.subplots()
        else:
            print args, type(args)
            fig, ax = args[0], args[1]
            x, y, z =self.extract_xyz()
            ax.plot(x, y)
        return fig, ax

    def analytical_singlepoint(self, velocity, direction=positions.CartesianPoint(1,0,0)):
        self.exact_solution = exact_translation(self.initial_position, velocity, direction)
        assert(self.exact_solution>0)

    def analytical_raypath(self, velocity, direction):
        pass


if __name__ == '__main__':

    
    print " test A " 
    A = Model()
    print A.Npoints
    A.initialisation(1., 10.,1.) 
    #plt.show()
    A.add_points(10.)
    print A.Npoints
    #plt.show()

    r, t, p = A.extract_rtp()
    plt.plot(r/max(r), '.')
    plt.plot(np.cos(np.pi/2-t*np.pi/180.), '.')
    #plt.show()
    plt.plot(p/max(p), '.')

    print "test B" 

    fig, ax = plt.subplots(1, 2)
    B = Translation()
    print "init"
    B.initialisation(0., 1221., 1e-6)
    print "init complete. ", B.Npoints, " points."
    B.tmax = 2.1 
    B.dt = 0.1
    B.vt =  1221.
    B.density_surface = 0.0010
    B.direction = np.array([1,0,0])
    x, y, z = B.extract_xyz()
    ax[0].plot(x, 'b.')
    ax[1].plot(x, 'b.')
    B.run(B.tmax, B.dt)
    
    x, y, z = B.extract_xyz()
    ax[0].plot(x, 'r.')
    print B.proxy
    #plt.show()


    
    r0, t0, t1, vt, omega = [.0,0.,0.], 0., 1., [1., 1., 1.], 0.
    rICB = 1221.
    print r0
    print integration_trajectory(r0, t0, t1, vt, omega)

    fig, ax = plt.subplots(1,2)
    ax[0].set_aspect('equal')

    point = 1221.*np.array([ -0.5, -0.5, 0.])
    t0, tmax, N = -2., 2., 30.
    vt = 1221.*np.array([1., 0.,     0. ])
    omega = 0.

    for t in np.linspace(t0, tmax, N):
        pos = integration_trajectory(point, t0, t, vt, omega)
        print "pos", pos
        r = trajectory_r(point, t0, t, vt, omega)
        ax[0].scatter(pos[0], pos[1])
        ax[1].scatter(t, r)

    solution = fsolve(lambda x : integration_trajectory(point, t0, x, vt, omega), -0.05)
    print solution
    solution = fsolve(lambda x : trajectory_r(point, t0, x, vt, omega)-rICB, -0.05)
    print solution
    x = np.linspace(-rICB, rICB , 100)
    ax[0].plot(x, np.sqrt(rICB**2-x**2), 'k')
    ax[0].plot(x, -np.sqrt(rICB**2-x**2), 'k')
    plt.show()
