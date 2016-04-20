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


class ModelGeodynamic():
    
    def __init__(self):
        self.rICB = 1221. #inner core radius in km.

    def set_tauIC(self, tau):
        self.tau_ic = tau

    def velocity(self, position, t):
        """ Velocity at the given position and given time. 
            
            Need to be implemented in derived classes. 
            Needed in cartesian coordinates.
            
            Example :
            velocity, in cartesian coordinates
    
            r is a cartesian position [x, y, z]
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
        solution = self.find_intersection(r0, t0, tau_ic)
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

    def find_intersection(self, r0, t0, x0):
        """ intersection between the trajectory and the radius, using fsolve method of scipy.optimize
            
            return the time corresponding to the intersection.
            x0 is the value used for starting the search (fsolve gives only one root)
            """
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

    def velocity(self, position, t):
        return self.vt

    def radius_ic(self, t):
        return self.rICB



if __name__ == '__main__':

    
    rICB = 1221.
    vt = 1221.*np.array([1., 0.,     0. ])
    t0 = 1.
    Model = PureTranslation(vt)
    Model.set_tauIC(t0)
    

    fig2, ax2 = plt.subplots(1,2)
    ax2[0].set_aspect('equal')

    point = 1221.*np.array([ 0.5, 0.5, 0.])
    point = [-1092.47368421, -449.842105263, 0.0]
    print point
    print point[0], point[1]
    
    tmax = 2.
    N = 10.
    for t in np.linspace(t0, tmax, N):
        pos = Model.integration_trajectory(t, point, t0)
        r = Model.trajectory_r(t, point, t0)
        ax2[0].scatter(pos[0], pos[1], c=t)
        ax2[1].scatter(t, r)
        print "time: ", t, "pos: ", pos

    solution = fsolve(lambda x : Model.integration_trajectory(x, point, t0), -0.05)
    print solution
    pos = Model.integration_trajectory(solution, point, t0)
    ax2[0].scatter(pos[0], pos[1], c='y')
    ax2[1].scatter(solution, np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2), c='y')
    solution = fsolve(lambda x : Model.trajectory_r(x, point, t0)-rICB, 2)
    print solution
    pos = Model.integration_trajectory(solution, point, t0)
    ax2[0].scatter(pos[0], pos[1], c='g')
    ax2[1].scatter(solution, np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2), c='g')
    u = np.linspace(-rICB, rICB , 100)
    ax2[0].plot(u, np.sqrt(rICB**2-u**2), 'k')
    ax2[0].plot(u, -np.sqrt(rICB**2-u**2), 'k')
    ax2[0].scatter(point[0], point[1], c='r')
    ax2[1].scatter(t0, np.sqrt(point[0]**2+point[1]**2+point[2]**2), c='r')
    
    print Model.find_intersection(point, t0, t0)
    print Model.find_time_beforex0(point, t0, t0)

    plt.show()
