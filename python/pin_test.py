import pinocchio
import numpy as np

urdf_filename = "./urdf/rm_75_6fb_description/urdf/RM75-6F.urdf"
model = pinocchio.buildModelFromUrdf(urdf_filename)
data = model.createData()
print(model)

# q = pinocchio.neutral(model)
q = np.array([0,1.58,0,0,0,0,0])

v = pinocchio.utils.zero(model.nv)
v[1] = -10
# print(v)
a = pinocchio.utils.zero(model.nv)

model.gravity[5] = 0
tau = pinocchio.rnea(model, data, q, v, a)
print("tau = ", tau)