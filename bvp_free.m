function [P, model, result] = bvp_free(C1, C2, Ri, Ro, rho, nu, varargin)
    %% Function for pressure in a contact-free bvp problem
    % constants C1, C2 = Material Moduli: Young's modulus and Bulk Modulus
    % P = pressure
    % rho = mass density of material
    % nu = poisson's ratio
    % varargin is for when ri and ro are given
    % ri = varargin{1}
    % ro = varargin{2};
    if length(varargin)==2
        ri = varargin{1};
        ro = varargin{2};
    elseif length(varargin)==1
        ri = varargin{1};
        ro = (Ro^3 + ri^3 - Ri^3)^(1/3);
    else
        disp('you must supply at least the internal radius, r_i');
    end

    gm = multisphere([Ri, Ro]);
    model = createpde('structural','static-solid');
    model.Geometry = gm;
% 
%     pdegplot(model,'CellLabels','on','FaceAlpha',0.4,'FaceLabels','on')
%     title('IAB with incompressibility constraints');

    % For each SoRo, specify the Young's modulus, Poisson ratio and mass
    % density
    structuralProperties(model, 'Cell', 1, 'YoungsModulus',C1,...
                'PoissonsRatio',nu,'MassDensity',rho);
    structuralProperties(model, 'Cell', 2, 'YoungsModulus',C2,...
                'PoissonsRatio',nu,'MassDensity',rho);
    % define the boundary conditions
    structuralBC(model,'Face',[1, 2],'Constraint','symmetric'); % symmetric isochoric

    %Specify the gravity load on the sphere.
    structuralBodyLoad(model, 'GravitationalAcceleration',[0;0;-9.8]);

    % relationship between radii in both configurations
    % ro = (Ro^3 + ri^3 - Ri^3)^(1/3);

    %specify contact-free pressure as a function handle (see thesis)
    fun = @(R, r) 2.* C1*((1./r)-((R.^6)./(r.^7))) + ...
               2.* C2*((r.^5)./(R.^6) - (R.^10)./(r.^11));

    P = integral(@(R)fun(R, ro), Ri, Ro);
    % apply the new pressure in Pascals
    structuralBoundaryLoad(model,'Face',2,'Pressure', P)

    % now create a mesh to view the deformation: see https://www.mathworks.com/help/pde/examples/deflection-analysis-of-a-bracket.html
    %sphere_thickness = 1e-1;
    start = tic;
    generateMesh(model); %, 'Hmax', sphere_thickness);
%     pdeplot3D(model)
%     title('Spherical IAB with Quadratic Tetrahedron Elements')
    fprintf('Time to mesh: %f seconds', toc(start))

    % calculate the solution
    result = solve(model);
end

