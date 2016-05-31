#!/usr/local/bin/python 
# Author: Marine Lasbleis


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq

## dummies function for testing
def function_f(x, *args):
    return np.sin(x+1.) 

def function_g(x):
    return np.sin(x) 

def function_diff(x, *args ):
    return function_f(x, *args)- function_g(x)


## solver functions

def zero_fsolve(f, x0, *args):
    """ find the zero of the function f using fsolve. 

    Please be careful that this method is not a good one when the derivative of the function f changes a lot
    x0 is the first point tested by fsolve. In the IC code, please use either 0 or age_ic."""
    return fsolve(lambda x : f(x), x0)

def zero_newton():
    pass #TODO

def zero_brentq(f, *args, **kwargs):
    """

    two possibilities: 
        - either the interval a,b is already set (and we just have to check the signs of f(a) and f(b).
        - or it is not, and a figure pops up to define the interval. 

        the program can also learn the value of the interval from other already calculated values. But need to be a little bit more careful (need to check the closests points, try to use this value, etc.)

    """
    if 'a' in kwargs and 'b' in kwargs:
        # simple one: the interval [a,b] is known
        a = kwargs["a"]
        b = kwargs["b"]
    else: 
        a, b = choose_interval_on_graph(f, 0, 2, *args) 
    
    if b>a: b,a = a,b

    iterations = 0
    while np.sign(f(a, *args)) == np.sign(f(b, *args)) and iterations <= 10 :
        if iterations == 0:
            print "Please choose other values for the interval so that the function changes sign between the two borns."
            a, b = choose_interval_on_graph(f, a, b, *args)
        elif iterations == 10:
            print "To solve this, you need to choose the interval such as the function changes sign. You have a last chance to find one: please type the two values of the interval."
            a = float(input("Enter first value: "))
            b = float(input("Enter second value: "))
        else:
            a, b = choose_interval_on_graph(f, a-(b-a), b+(b-a), *args)
        iterations += 1    

    return brentq(f, a, b, args = args)

def choose_interval_on_graph(f, a, b, *args):
    """ allow the user to click on graph to choose an interval."""
    global coords
    coords = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_val = np.linspace(a, b, 100)
    f_val = np.zeros_like(x_val)
    for i, x in enumerate(x_val): f_val[i] = f(x, *args)
    ax.plot(x_val, f_val) 
    ax.plot(x_val, x_val*0.)
    ax.set_title("Please click on two points to define the interval.")
    cid = fig.canvas.mpl_connect('button_press_event', onclick_two)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    return coords[0], coords[1] 



# 
# def intersection(h, g, x0, *args, **kwargs):
#     """ find the intersection of function h and g, starting looking at x0"""
#     if kwargs:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         print "enter kwargs"
#         print kwargs
#         if ("min" in kwargs) and ("max" in kwargs):
#             x_val = np.linspace(kwargs["min"], kwargs["max"], 100)
#             ax.plot(x_val, h(x_val, *args))
#             ax.plot(x_val, g(x_val))
#             func = lambda x : h(x, *args) - g(x)
#             ax.plot(x_val, func(x_val) )
#             cid = fig.canvas.mpl_connect('button_press_event', onclick)
#             plt.show()
#             fig.canvas.mpl_disconnect(cid)
#             print position_x
#             solution = fsolve(func, position_x) 
#             print solution
#             # TODO use the brentq formulat to obtain the roots!
#         else: print "problem. Please give values for min and max"
#             #do a figure check (plot the figure of both curves + ask the person to select where he wants to do the search)
#     else: solution = zero_brentq(function_diff, a=5)#, b=3) 
#     # fsolve(lambda x : h(x, *args) - g(x), x0)
#     return  solution 

def onclick(event):
    global position_x
    print('button=%d, xdata=%f, ydata=%f' %
          (event.button, event.xdata, event.ydata))
    position_x = event.xdata
    plt.close()

def onclick_two(event):
    global coords
    print('button=%d, xdata=%f, ydata=%f' %
          (event.button, event.xdata, event.ydata))
    ix = event.xdata
    coords.append(ix)
    if len(coords)==2:
        plt.close()


if __name__ == "__main__":

    c = intersection(function_f, function_g, 20)
    print   c, function_f(c), function_g(c), function_diff(c) 

