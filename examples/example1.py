

import sys
# in this case local import
sys.path.append("../")
import mag2dpoly as mag

import numpy as np


# induced magnetization
Jind = mag.MagnetizVector(mod=4.9,Ideg=90.0,Ddeg=45.0)
# remanent magnetization
Jrem = mag.MagnetizVector(mod=3.1,Ideg=45.0,Ddeg=0.0)

# angle with the North axis
northxax = 90.0

# number of observation
Nobs = 101
xzobs = np.transpose(np.vstack(( np.linspace(0.0,100.0,Nobs), -1.0*np.ones(Nobs))))

# vertices of the poligonal bodies
vertices = np.array([ [35.0, 50.0],
                      [65.0, 50.0],
                      [80.0, 35.0],
                      [65.0, 20.0],
                      [35.0, 20.0],
                      [20.0, 35.0] ])

# indices of vertices for the body
nbod = 1
bodyindices = np.empty(shape=(nbod,), dtype=np.object)
inds = range(6)
bodyindices[0] = np.array(inds)

# construct the poligonal body object
pbody = mag.MagPolyBodies2D(bodyindices,vertices)


# type of forward algorithm
forwardtype = "talwani"

# compute total field
# make Jind and Jrem arrays of objects (as many as there are bodies)
Jindv = np.array([Jind]) # we have one single body in this case
Jremv = np.array([Jrem]) # we have one single body in this case
tmag = mag.tmagpolybodies2Dgen(xzobs,Jindv,Jremv,northxax,pbody,forwardtype)


## plot
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(211)
plt.title("Magnetic anomaly")
plt.plot(xzobs[:,0],tmag,"o-")
plt.subplot(212)
plt.title("Polygonal body")
x = np.append(pbody.bo[0].ver1[:,0],pbody.bo[0].ver1[0,0])
y = np.append(pbody.bo[0].ver1[:,1],pbody.bo[0].ver1[0,1])
plt.plot(x,y,"o-")
plt.show()
