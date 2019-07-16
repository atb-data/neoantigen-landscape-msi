import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def cmap_map(function, cmap):
  """ Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
  This routine will break any discontinuous points in a colormap.
  """
  cdict = cmap._segmentdata
  step_dict = {}
  # First get the list of points where the segments start or end
  for key in ('red', 'green', 'blue'):
    step_dict[key] = list(map(lambda x: x[0], cdict[key]))
  step_list = sum(step_dict.values(), [])
  step_list = np.array(list(set(step_list)))
  # Then compute the LUT, and apply the function to the LUT
  reduced_cmap = lambda step : np.array(cmap(step)[0:3])
  old_LUT = np.array(list(map(reduced_cmap, step_list)))
  new_LUT = np.array(list(map(function, old_LUT)))
  # Now try to make a minimal segment definition of the new LUT
  cdict = {}
  for i, key in enumerate(['red','green','blue']):
    this_cdict = {}
    for j, step in enumerate(step_list):
      if step in step_dict[key]:
        this_cdict[step] = new_LUT[j, i]
      elif new_LUT[j,i] != old_LUT[j, i]:
        this_cdict[step] = new_LUT[j, i]
    colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
    colorvector.sort()
    cdict[key] = colorvector

  return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
  new_cmap = colors.LinearSegmentedColormap.from_list(
    'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
    cmap(np.linspace(minval, maxval, n))
  )
  return new_cmap

def desaturate(cmap):
  def desaturate_impl(vec):
    return (np.array([vec.sum() / 3.0] * 3) * 5.0 / 8.0 + vec * 3.0 / 8.0)
  return cmap_map(desaturate_impl, cmap)

def cmap_discretize(cmap, N):
  """Return a discrete colormap from the continuous colormap cmap.

  Args:
    cmap: colormap instance, eg. cm.jet.
    N: number of colors.
  """
  colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
  colors_rgba = cmap(colors_i)
  indices = np.linspace(0, 1., N+1)
  cdict = {}
  for ki, key in enumerate(('red','green','blue')):
    cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1)]
  # Return colormap object.
  return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
