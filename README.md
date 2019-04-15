Some useful formulae in screw theory for robotics

Detailed README coming soon, posr-thesis defense

+ [rodrigues](/rodrigues.m): Given a vector, w, and an angle, theta, computes the exponential map of an angle using Rodrigues' formula. This is akin to the rigid body transformation of a body about an axis w, through a point q, by an angle theta.

+ [skewsem](skewsem.m): Computes the antisymmetric or skew symmetric matrix from a 3D vector w = [wx, wy, wz].

+ [twist](/twist.m): Given an infinitesimal linear generator v and an infinitesimal angular generator w, computes the twist vector \in R^6 as [(-w x v) w]

+ [twistmap](/twistmap.m): Given a twist, \zeta, and an angle of rotation \theta, computes the exponential map of the twist exp(\hat{zeta}, theta), where $\hat{\zeta}$ is the isomorphism of the twist linear algebra in R^6

+ [curvature_form](/curvature_form.m): Computes the relative curvature at the point of contact of two objects. This is typically useful in multifingered kinematics or IAB manipulations such as I used in my PhD thesis. This implementation is from Spivak's differential geometry book, 1979.

+[psi_dot](/psi_dot.m): computes the relative orientation between the tangent planes of two coordinate charts in R^2. The rotation angle through an angle \psi aligns the manipulating object's x-axis to the x-axis of the manipulated object. \psi is so chosen such that a rotation about the z-axis by \psi aligns the xy-axes of both objects (see Montana, 1988).

+ [contact_kinematics](/contact_kinematics.m): A simple illustration of the contact dynamics between two rigid, or semi rigid objects in R^3. The examples assume that the objects are spheres with single coordinate charts as defined in [Murray, 1990].
