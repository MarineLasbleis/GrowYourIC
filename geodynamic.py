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
            time[i] = evaluate_singlepoint(point, method)[0]
        elif dataset.method == "raypath":
            N = dataset.NpointsRaypath
            dataset.data_points[i].straigth_in_out(N)
            raypath = ray.points
            total_proxy = 0.
            for j, point in enumerate(raypath):
                _proxy = evaluate_singlepoint(point, method)[0]
                total_proxy += _proxy
            time[i] = total_proxy / float(N)
    return time

def evaluate_singlepoint(point, method):
    """ evaluate the proxy on a single positions.Point instance, using the choosen method."""
    x, y, z = point.x, point.y, point.z
    time = method.find_time_beforex0([x, y, z], method.tau_ic, method.tau_ic)
    return method.tau_ic-time

def trajectory_single_point(point, method, t0, t1, num_t):
    """ return the trajectory of a point (a positions.Point instance) between the times t0 and t1, knowing that it was at the position.Point at t0, given nt times steps. 
    """
    time = np.linspace(t0, t1, num_t)
    x, y, z = np.zeros(num_t), np.zeros(num_t), np.zeros(num_t)
    x[0], y[0], z[0] = point.x, point.y, point.z
    for i, t in enumerate(time):
        point = method.integration_trajectory(t, [x[0], y[0], z[0]], t0)
        x[i], y[i], z[i] = point[0], point[1], point[2]

    return x, y, z 


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

    def set_tauIC(self, tau):
        self.tau_ic = tau

    def set_exponent_growth(self, alpha):
        self.exponent_growth = alpha
    
    def set_rICB(self, RIC):
        self.rICB = RIC#value by default is 1221, but can be changed if necessary. 

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


    def find_time_beforex0(self, r0, t0, tau_ic, dx=1.):
        """ find the first intersection between the trajectory and the radius of the IC
            
            tau_ic is the age of inner core, and only times of intersection smaller than x0 are output.
            """
        solution = self.find_intersection(r0, t0, 0.) #tau_ic)
        if tau_ic != 0:
            dx = max(tau_ic/300., min(dx, tau_ic/20.))
            # just to check if the dx has been chosen OK. Assuming tau_ic is the age of the IC, dx should be about tau_ic/50. It should not be bigger than tau_ic/50 or smaller than x0/300. (rule of thumb)
        N =0
        while solution > tau_ic  and N <= 300:
            N += 1
            solution = self.find_intersection(r0, t0, tau_ic-N*dx)
            if N == 300:
                print "No intersection found."
        return solution

    def find_intersection(self, r0, t0, x0, test=False):
        """ intersection between the trajectory and the radius, using fsolve method of scipy.optimize
        
        return the time corresponding to the intersection.
        x0 is the value used for starting the search (fsolve gives only one root)
        """
        if test:    
            time = np.linspace(0., self.tau_ic, 20)
            trajectory = np.zeros_like(time)
            radius = np.zeros_like(time)
            for i, t in enumerate(time):
                trajectory[i] = self.trajectory_r(t, r0, t0)
                radius[i] = self.radius_ic(t)
            solution = fsolve(lambda x : self.trajectory_r(x, r0, t0)-self.radius_ic(x), x0)     
            plt.plot(time, trajectory, time, radius)
            plt.scatter(solution, self.radius_ic(solution))
            plt.show()
        return fsolve(lambda x : self.trajectory_r(x, r0, t0)-self.radius_ic(x), x0)

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


class PureTranslation(ModelGeodynamic):
    
    def __init__(self, vt):
        ModelGeodynamic.__init__(self)
        self.name = "Translation"
        self.vt = vt # a np array (1,3)
        assert (len(vt)==3), "The velocity need to be specified as [Vx, Vy, Vz]"

    def velocity(self, t, position):
        return self.vt

    def radius_ic(self, t):
        return self.rICB

class TranslationRotation(ModelGeodynamic):

    def __init__(self, vt, omega):
        ModelGeodynamic.__init__(self)
        self.name = "TranslationRotation"
        self.vt = vt # a np array (1,3)
        self.omega = omega
    
    def velocity(self, t, r):
        """ velocity at the point position

        position is a np.array [x, y, z] 
        """
        
        return self.vt + np.array([-self.omega * r[1], self.omega * r[0], 0.] )
   
    def radius_ic(self, t):
        return self.rICB


class PureRotation(ModelGeodynamic):

    def __init__(self, omega):
        ModelGeodynamic.__init__(self)
        self.name = "Rotation"
        self.omega = omega
    
    def velocity(self, t, r):
        """ velocity at the point position

        position is a np.array [x, y, z] 
        """
        
        return np.array([-self.omega * r[1], self.omega * r[0], 0.] )
   
    def radius_ic(self, t):
        return self.rICB



class PureGrowth(ModelGeodynamic):
    def __init__(self):
        ModelGeodynamic.__init__(self)
        self.name = "PureGrowth"
        self.exponent_growth = None 
    
    def velocity(self, t, r):
        """ velocity at the point position
    
        position is a np.array [x, y, z]
        """
        return np.array([0.,0.,0.]) 
  
    def radius_ic(self, t):
        return self.rICB*(t/self.tau_ic)**self.exponent_growth



class TranslationGrowth(ModelGeodynamic):

    def __init__(self, vt):
        ModelGeodynamic.__init__(self)
        self.name = "Translation and Growth"
        self.exponent_growth = None 
        self.vt = vt # a np array (1,3)


    def velocity(self, t, r):
        """ velocity at the point position
    
        position is a np.array [x, y, z]
        """
        return self.vt 
  
    def radius_ic(self, t):
        return self.rICB*(t/self.tau_ic)**self.exponent_growth




if __name__ == '__main__':

    vt = [1.,0.,0.]
    omega = np.pi 
    Method = TranslationGrowth(vt)
    Method = TranslationRotation(vt,omega)
    Method.set_tauIC(1.)
    Method.set_exponent_growth(0.5)
    Method.set_rICB(1.)

    point = positions.CartesianPoint(0.5, -0.5, 0.)
    traj_x, traj_y, traj_z = trajectory_single_point(point, Method,  1., 0., 10 )
    print traj_x, traj_y, traj_z
    plt.plot(traj_x, traj_y, label=Method.name)

    Method = PureTranslation(vt)
    traj_x, traj_y, traj_z = trajectory_single_point(point, Method,  1., 0., 10 )
    plt.plot(traj_x, traj_y, label=Method.name)

    Method = PureGrowth()
    traj_x, traj_y, traj_z = trajectory_single_point(point, Method,  1., 0., 10 )
    plt.plot(traj_x, traj_y,'o-',  label=Method.name)


    Method = PureRotation(omega)
    traj_x, traj_y, traj_z = trajectory_single_point(point, Method,  1., 0., 10 )
    plt.plot(traj_x, traj_y,'-',  label=Method.name)


    phi = np.linspace(0, 2*np.pi)
    x = np.cos(phi)
    y = np.sin(phi)
    plt.plot(x,y)

    plt.axis('equal')
    plt.legend()
    plt.show()
