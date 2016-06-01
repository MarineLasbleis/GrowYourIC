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
import intersection 

RICB = 1221.

def evaluate_proxy(dataset, method):
    """ evaluate the value of the proxy on all the points of the data set, using the choosen geodynamical method
        
        dataset : a data.SeismicData object
        method : a geodynamic.ModelGeodynamic object
        """
    # TO DO : choose evaluate on raypath or on BT point
    print "==="
    print "== Evaluate value of proxy for all points of the data set "
    print "= Geodynamic model is ", method.name
    print "= Data set is ", dataset.name
    print "= Proxy is evaluate for ", dataset.method
    print "= Number of points to examine: ", dataset.size 

    time = np.empty_like(dataset.data_points)
    for i, ray in enumerate(dataset.data_points):
        if dataset.method == "bt_point":
            point = ray.bottom_turning_point
            time[i] = method.proxy_singlepoint(point)
        elif dataset.method == "raypath":
            N = dataset.NpointsRaypath
            dataset.data_points[i].straigth_in_out(N)
            raypath = ray.points
            total_proxy = 0.
            for j, point in enumerate(raypath):
                _proxy = method.proxy_singlepoint(point)
                total_proxy += _proxy
            time[i] = total_proxy / float(N)
    return time 



def exact_translation(point, velocity, direction=positions.CartesianPoint(1,0,0)):
    x_0, y_0, z_0 = point.x, point.y, point.z
    mean_direction = np.sqrt(direction.x**2+ direction.y**2+ direction.z**2)
    a, b, c = direction.x/mean_direction, direction.y/mean_direction, direction.z/mean_direction

    solution_1 = x_0*a+y_0*b+z_0*c + np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**2))
    solution_2 = x_0*a+y_0*b+z_0*c - np.sqrt((x_0*a+y_0*b+z_0*c)**2-(x_0**2+y_0**2+z_0**2-RICB**    2))
    ## TO DO : verify that we can remove solution_2 ? Is solution_1 always max?
    return max(solution_1, solution_2)


class ModelGeodynamic():
    
    def __init__(self):
        self.rICB = 1221. #inner core radius in km.

    def set_parameters(self, dict_param):
        """ add any parameter of the form {'param':value} as self.param = value """
        for k, v in dict_param.items():
            setattr(self, k, v)

    def set_tauIC(self, tau):
        self.tau_ic = tau

    def set_exponent_growth(self, alpha):
        self.exponent_growth = alpha
    
    def set_rICB(self, RIC):
        self.rICB = RIC#value by default is 1221, but can be changed if necessary. 

    def set_vt(self, vt):
        while True:
            try:
                assert(len(vt)==3)
            except AssertionError:
                print "Translation velocity needs to have 3 components. Please enter the 3 cartesians components of the velocity: "
                vt = map(float, raw_input("Translation velocity: ").split())
                continue #come back to check if valid answer
            else:
                break #no error raise!
        self.vt = vt

    def set_rotation(self, omega):
        self.omega = omega


    def velocity(self, t, position):
        """ Velocity at the given position and given time. 
            
            Need to be implemented in derived classes. 
            Needed in cartesian coordinates.
            
            Example :
            velocity, in cartesian coordinates
    
            position is a cartesian position [x, y, z]
            vt is the translation velocity (numpy array (1,3))
            rotation velocity is $\omega \vec{e}_z \times \vec{r}$.
            
            return vt + np.array([-omega * r[1], omega * r[0], 0.] )
            """
        raise NotImplementedError("need to implement velocity() in derived class!")

    def radius_ic(self, t):
        """ radius of the inner core with time. 
            
            Need to be implemented in derived classes.
            """
        raise NotImplementedError("need to implement velocity() in derived class!")


    def find_time_beforex0(self, point, t0, t1):
        """ find the intersection between the trajectory and the radius of the IC
        if needed, can be re defined in derived class!

        point : [x, y, z]
        """
        return intersection.zero_brentq(self.distance_to_radius, point, t0, a=0., b=t1)


    def proxy_singlepoint(self, point):
        """ evaluate the proxy on a single positions.Point instance."""
        ## TODO proxy here is age only. Please change this if needed.
        if point.r< self.rICB:
            x, y, z = point.x, point.y, point.z
            time = self.find_time_beforex0([x, y, z], self.tau_ic, self.tau_ic)
        else:
            x, y, z = point.x, point.y, point.z
            time = self.find_time_beforex0([x, y, z], self.tau_ic, self.tau_ic*1.01)
        
        return self.tau_ic-time

    def distance_to_radius(self, t, r0, t0):
        return self.trajectory_r(t, r0, t0)-self.radius_ic(t)

    def trajectory_r(self, t, r0, t0):
        """ for a point at position r0 at time t0, return the radial component of the position of the point at time t.
            
            """
        trajectory = self.integration_trajectory(t, r0, t0)
        #r, t, p = positions.from_cartesian_to_seismo(trajectory[0], trajectory[1], trajectory[2])
        r = trajectory[0]**2+trajectory[1]**2+trajectory[2]**2
        return np.sqrt(r)

    def integration_trajectory(self, t1, r0, t0):
        """ integration of the equation dr(t)/dt = v(r,t)
    
        r0: initial position
        t0: initial time
        t1: tmax of the integration
            """
        r = ode(self.velocity).set_integrator('dopri5')
        r.set_initial_value(r0, t0) # .set_f_params() if the function has any parameters
        return np.real(r.integrate(r.t+(t1-t0)))

    def trajectory_single_point(self, point, t0, t1, num_t):
        """ return the trajectory of a point (a positions.Point instance) between the times t0 and t1, knowing that it was at the position.Point at t0, given nt times steps. 
        """
        time = np.linspace(t0, t1, num_t)
        x, y, z = np.zeros(num_t), np.zeros(num_t), np.zeros(num_t)
        x[0], y[0], z[0] = point.x, point.y, point.z
        for i, t in enumerate(time):
            point = self.integration_trajectory(t, [x[0], y[0], z[0]], t0)
            x[i], y[i], z[i] = point[0], point[1], point[2]
        return x, y, z 
    
    def plot_equatorial(self, t0, t1, Nt = 200,  N=40):
        # Plot the inner core boundary
        phi = np.linspace(0., 2*np.pi, N)
        x, y = self.rICB*np.cos(phi), self.rICB*np.sin(phi)
        plt.plot(x, y, 'b')
        for i, phi in enumerate(phi):
            trajx, trajy, trajz = self.trajectory_single_point(positions.CartesianPoint(x[i], y[i], 0.), t0, t1, Nt)
            trajectory_r = np.sqrt(trajx**2. + trajy**2 + trajz**2)
            mx = np.ma.masked_array(trajx, mask = trajectory_r >= self.rICB)
            my = np.ma.masked_array(trajy, mask = trajectory_r >= self.rICB)
            plt.plot(mx, my)
            plt.plot(mx[::10], my[::10], '+b')
            velocity = self.velocity(self.tau_ic, [trajx, trajy, trajz])
            plt.quiver(mx[::10], my[::10], np.ma.masked_array(velocity[0], mask = trajectory_r >= self.rICB)[::10], np.ma.masked_array(velocity[1], mask = trajectory_r >= self.rICB)[::10], units='width')
        plt.axis("equal")
        plt.xlim([-1,1])
        plt.ylim([-1,1])
        #plt.show()

    def translation_velocity(self):
        try:
            self.vt
            assert(len(self.vt)==3)
            return self.vt
        except (NameError, AttributeError):
            print "translation velocity has not been given. Please enter the three components of the velocity: (be careful, it has to be dimensionless)"
            value = map(float, raw_input("Translation velocity: ").split())
            self.set_vt(value)

    def rotation_velocity(self, r):
        try:
            self.omega
        except (NameError, AttributeError):
            print "Rotation velocity has not been defined. Please enter the value for omega: "
            value = float(input("Rotation rate: "))
            self.set_rotation(value)
        return np.array([-self.omega * r[1], self.omega * r[0], 0.] )    

    def growth_ic(self, t):
        try:
            self.rICB
        except (AttributeError, NameError):
            print "The value of rICB has not been provided. Please enter it now: "
            value = float(input("rICB: "))
            self.set_rICB(value)
        try:
            assert(self.tau_ic != None)
        except (AttributeError, NameError, AssertionError):
            print "The value of tau_ic has not been provided. Please enter it now: "
            value = float(input("tau_ic: "))
            self.set_tauIC(value)
        try:
            assert(self.exponent_growth != None)
        except (AttributeError, NameError, AssertionError):
            print "The value of exponent_growth has not been provided. Please enter it now: "
            value = float(input("exponent growth: "))
            self.set_exponent_growth(value)
        return self.rICB*(t/self.tau_ic)**self.exponent_growth

class PureTranslation(ModelGeodynamic):
    
    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "Translation"

    def velocity(self, t, position):
      #  position = np.array(position)
      #  if position.ndim == 1: 
      #      ncols= 1
      #  else: 
      #      nlines, ncols = position.shape
      #  return np.repeat([self.vt], ncols, axis=0)
        return self.translation_velocity() 

    def radius_ic(self, t):
        return self.rICB
    



class TranslationRotation(ModelGeodynamic):

    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "TranslationRotation"
    
    def velocity(self, t, r):
        """ velocity at the point position

        position is a np.array [x, y, z] 
        """
        return self.translation_velocity() + self.rotation_velocity(r)
   
    def radius_ic(self, t):
        return self.rICB

class PureRotation(ModelGeodynamic):

    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "Rotation"
    
    def velocity(self, t, r):
        """ velocity at the point position

        position is a np.array [x, y, z] 
        """
        
        return self.rotation_velocity() 
   
    def radius_ic(self, t):
        return self.rICB



class PureGrowth(ModelGeodynamic):

    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "PureGrowth"
    
    def velocity(self, t, r):
        """ velocity at the point position
    
        position is a np.array [x, y, z]
        """
        return np.array([0.,0.,0.]) 
  
    def radius_ic(self, t):
        return self.growth_ic(t) 



class TranslationGrowth(ModelGeodynamic):

    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "Translation and Growth"

    def velocity(self, t, r):
        """ velocity at the point position
    
        position is a np.array [x, y, z]
        """
        return self.translation_velocity() 
  
    def radius_ic(self, t):
        return self.growth_ic(t)


class TranslationGrowthRotation(ModelGeodynamic):

    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "Translation, Rotation and Growth"


    def velocity(self, t, r):
        """ velocity at the point position
    
        position is a np.array [x, y, z]
        """
        return self.translation_velocity()+ self.rotation_velocity(r)
  
    def radius_ic(self, t):
        return self.growth_ic(t)





if __name__ == '__main__':

    vt = [2.,0.,0.]
    omega = 0.5*np.pi 

    points = [positions.CartesianPoint(-0.2, -0.85, 0.)]
    N = 20

    for point in points:


        Method = TranslationRotation(vt,omega)
        Method.set_tauIC(1.)
        Method.set_exponent_growth(0.5)
        Method.set_rICB(1.)
        traj_x, traj_y, traj_z = Method.trajectory_single_point(point,  1., 0., N )
        plt.plot(traj_x, traj_y, label=Method.name)

        time = Method.find_intersection([point.x, point.y, point.z], 1., 0.)
        pos = Method.integration_trajectory(time, [point.x, point.y, point.z], 1.)
        plt.plot(pos[0], pos[1], 'bo')
    
        Method = PureTranslation(vt)
        Method.set_tauIC(1.)
        Method.set_exponent_growth(0.5)
        Method.set_rICB(1.)
 
        traj_x, traj_y, traj_z = Method.trajectory_single_point(point, 1., 0., N )
        plt.plot(traj_x, traj_y, label=Method.name)
        time = Method.find_intersection([point.x, point.y, point.z], 1., 0.)
        pos = Method.integration_trajectory(time, [point.x, point.y, point.z], 1.)
        plt.plot(pos[0], pos[1], 'bo')

        Method = PureGrowth()
        Method.set_tauIC(1.)
        Method.set_exponent_growth(0.5)
        Method.set_rICB(1.)
        traj_x, traj_y, traj_z = Method.trajectory_single_point(point,  1., 0., N )
        plt.plot(traj_x, traj_y,'o-',  label=Method.name)
        time = Method.find_intersection([point.x, point.y, point.z], 1., 3.)
        pos = Method.integration_trajectory(time, [point.x, point.y, point.z], 1.)
        plt.plot(pos[0], pos[1], 'bo')



        Method = PureRotation(omega)
        Method.set_tauIC(1.)
        Method.set_exponent_growth(0.5)
        Method.set_rICB(1.)
        traj_x, traj_y, traj_z = Method.trajectory_single_point(point,  1., 0., N )
        plt.plot(traj_x, traj_y,'-',  label=Method.name)


    phi = np.linspace(0, 2*np.pi, 50)
    x = np.cos(phi)
    y = np.sin(phi)
    plt.plot(x,y)

    plt.axis('equal')
    plt.legend()
    plt.show()

    positions_x = np.linspace(-1., 1., 10)
    positions_y = np.linspace(-1., 1, 10)
    points = []
    for x in positions_x:
        for y in positions_y:
            points.append(positions.CartesianPoint(x, y, 0.)) 
    
    N = 20

    for point in points:


        Method = TranslationRotation(vt,omega)
        Method.set_tauIC(1.)
        Method.set_exponent_growth(0.5)
        Method.set_rICB(1.)
        traj_x, traj_y, traj_z = Method.trajectory_single_point(point,  1., 0., N )
        plt.plot(traj_x, traj_y, label=Method.name)



    phi = np.linspace(0, 2*np.pi, 50)
    x = np.cos(phi)
    y = np.sin(phi)
    plt.plot(x,y)

    plt.axis('equal')
    plt.show()

