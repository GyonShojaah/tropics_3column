import os, itertools, cPickle, hashlib
import scipy.interpolate
import numpy as np

##############################
## Interpolation
##############################
# Scipy has many functions to do 2d interpolation.  The all have
# different calling conventions.  Some of them aren't very good, some
# aren't robust, some don't respect all of the arguments you give
# them.  I want to feed three 2D arrays of x, y, and z values.  I want
# to optionally do the interpolation in linear or log space for each
# of the variables.  For the resulting function, I want the inputs and
# outputs to be in linear space regardless of how it does ithe
# interpolation.  Finally, I want it to apply pointwise to its
# arguments for scalars, 1d, and 2d arrays, rather than the wacky
# things that some of the scipy functions try to do.  This functions
# provide a uniform interface to the interpolation routines that I
# have found to be reasonably robust.

def interp_1d(xx,yy,logx=True, logy=True, kind='linear', bounds_error=True, fill_value=np.nan):
    xx = np.log(xx) if logx else xx
    yy = np.log(yy) if logy else yy
    internal = scipy.interpolate.interp1d(xx, yy, kind=kind, bounds_error=bounds_error, fill_value=fill_value)
    def interpolator(uu):
        uu = np.log(uu) if logx else uu
        ww = internal(uu)
        result = np.exp(ww) if logy else ww
        return result
    return interpolator


def interp_1d_boundary(xx,yy,logx=True, logy=True, order=1, bbox=[None, None]):
    xx = np.log(xx) if logx else xx
    yy = np.log(yy) if logy else yy
    internal = scipy.interpolate.InterpolatedUnivariateSpline(xx, yy, k=order)
    def interpolator(uu):
        uu = np.log(uu) if logx else uu
        ww = internal(uu)
        result = np.exp(ww) if logy else ww
        return result
    return interpolator
    

def interp_rect_spline(xx,yy,Z,logx=True, logy=True, logz=True, 
                       kx=3, ky=3, **kw):
    """interplation using splines on a rectangular grid."""
    # this one is different:  needs rectangular grid, so enforce this in the arguments.
    xx = np.log(xx) if logx else xx
    yy = np.log(yy) if logy else yy
    Z = np.log(Z) if logz else Z
    xx,yy,Z = [np.asarray(aa) for aa in xx,yy,Z]
    internal = scipy.interpolate.RectBivariateSpline(xx, yy, Z, 
                                                     kx=kx, ky=ky, **kw)
    def interpolator(uu,vv):
        uu,vv = np.asarray(uu), np.asarray(vv)
        if uu.shape != vv.shape: raise RuntimeError
        # interpolating function operates  on the cross product of
        # its  arguments, and  I  want it  to  operate element  by
        # elemnt.   This complicates  things...  Handle  the three
        # common cases: scalars, arrays, and 2d arrays
        if len(uu.shape) == 0:
            uu = np.log(uu) if logx else uu
            vv = np.log(vv) if logy else vv            
            zz = internal(uu,vv)[0,0]
            result = np.exp(zz) if logz else zz
        elif len(uu.shape) == 1: 
            uu = np.log(uu) if logx else uu
            vv = np.log(vv) if logy else vv            
            zz = np.array([internal(the_x, the_y)[0,0]
                           for the_x, the_y in zip(uu,vv)])
            result = np.exp(zz) if logz else zz
        elif len(uu.shape) == 2: 
            # now we're getting tricky
            result = interpolator(uu.ravel(), vv.ravel()).reshape(uu.shape)
        else: 
            raise RuntimeError 
        return result
    return interpolator

# Default 2d interpolation uses rect_spline.  
# interp_2d = interp_rect_spline

def interp_rbf(X,Y,Z,logx=True, logy=True, logz=True, **kw):
    """Unstructured interplation using radial basis functions"""
    # X,Y,Z should be physical values.  logx, etc, say whether you
    # want to do the interpolation in log space or not.
    X = np.log(X) if logx else X
    Y = np.log(Y) if logy else Y
    Z = np.log(Z) if logz else Z
    X,Y,Z = [np.asarray(aa) for aa in X,Y,Z]
    internal = scipy.interpolate.Rbf(X.ravel(), Y.ravel(), Z.ravel(), **kw)
    def interpolator(xx, yy):
        xx = np.log(xx) if logx else xx
        yy = np.log(yy) if logy else yy            
        zz = internal(xx, yy)
        return np.exp(zz) if logz else zz
    return interpolator


def interp_griddata(X,Y,Z,logx=True, logy=True, logz=True, **kw):
    """Unstructured interplation using nearest neighbors for now."""
    X,Y,Z = [np.asarray(aa) for aa in X,Y,Z]
    X = np.log(X) if logx else X
    Y = np.log(Y) if logy else Y
    Z = np.log(Z) if logz else Z
    coords = grid_to_points(X,Y)
    values = grid_to_points(Z)
    def interpolator(xx,yy):
        xx,yy = np.asarray(xx), np.asarray(yy)
        if xx.shape != yy.shape: raise RuntimeError
        xx = np.log(xx) if logx else xx
        yy = np.log(yy) if logy else yy            
        desired_coords = grid_to_points(xx,yy)
        # method=nearest is the only way I can get this to give sensible results.
        zz = scipy.interpolate.griddata(coords, values, desired_coords, 
                                        method='nearest', **kw)
        zz = zz.reshape(xx.shape)
        return np.exp(zz) if logz else zz
    return interpolator

def interp_griddata_3Dgrid(coords, values, points, logx=True, logy=True, logz=True, logw=True, **kw):
    """Unstructured interplation using nearest neighbors for now."""
    X,Y,Z,W = [np.asarray(aa) for aa in X,Y,Z,W]
    X = np.log(X) if logx else X
    Y = np.log(Y) if logy else Y
    Z = np.log(Z) if logz else Z
    W = np.log(Z) if logz else W
    coords = grid_to_points(X,Y)
    values = grid_to_points(Z)
    def interpolator(xx,yy):
        xx,yy = np.asarray(xx), np.asarray(yy)
        if xx.shape != yy.shape: raise RuntimeError
        xx = np.log(xx) if logx else xx
        yy = np.log(yy) if logy else yy            
        desired_coords = grid_to_points(xx,yy)
        # method=nearest is the only way I can get this to give sensible results.
        zz = scipy.interpolate.griddata(coords, values, desired_coords, 
                                        method='nearest', **kw)
        zz = zz.reshape(xx.shape)
        return np.exp(zz) if logz else zz
    return interpolator


def interp_spline(X,Y,Z,logx=True, logy=True, logz=True, **kw):
    """Unstructured interplation using splines."""
    X = np.log(X) if logx else X
    Y = np.log(Y) if logy else Y
    Z = np.log(Z) if logz else Z
    X, Y, Z = [np.asarray(aa) for aa in X,Y,Z]
    internal = scipy.interpolate.SmoothBivariateSpline(
        X.ravel(), Y.ravel(), Z.ravel(), **kw)
    def interpolator(xx,yy):
        xx,yy = np.asarray(xx), np.asarray(yy)
        if xx.shape != yy.shape: raise RuntimeError
        # interpolating function operates  on the cross product of
        # its  arguments, and  I  want it  to  operate element  by
        # elemnt.   This complicates  things...  Handle  the three
        # common cases: scalars, arrays, and 2d arrays
        if len(xx.shape) == 0:
            xx = np.log(xx) if logx else xx
            yy = np.log(yy) if logy else yy            
            zz = internal(xx,yy)[0,0]
            result = np.exp(zz) if logz else zz
        elif len(xx.shape) == 1: 
            xx = np.log(xx) if logx else xx
            yy = np.log(yy) if logy else yy            
            zz = np.array([internal(the_x, the_y)[0,0]
                           for the_x, the_y in zip(xx,yy)])
            result = np.exp(zz) if logz else zz
        elif len(xx.shape) == 2: 
            # now we're getting tricky
            result = interpolator(xx.ravel(), yy.ravel()).reshape(xx.shape)
        else: 
            raise RuntimeError 
        return result
    return interpolator
