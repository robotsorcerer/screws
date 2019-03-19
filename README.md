Some useful formulae in screw theory for robotics

Detailed README coming soon, posr-thesis defense

+ [rodrigues](/rodrigues.m): Given a vector, w, and an angle, theta, computes the exponential map of an angle. This is akin to the rigid body transformation of a body about an axis w, through a point q, by an angle theta

+ [skewsem](skewsem.m): Computes the antisymmetric or skew symmetric matrix from a 3D vector w = [wx, wy, wz]

+ [twist](/twist.m): Given an infinitesimal linear generator v and an infinitesimal angular generator w, computes the twist vector \in R^6 as [(-w x v) w]

+ [twistmap](/twistmap.m): Given a twist, \zeta, and an angle of rotation \theta, computes the exponential map of the twist exp(\hat{zeta}, theta), where $\hat{\zeta}$ is the isomorphism of the twist linear algebra in R^6
