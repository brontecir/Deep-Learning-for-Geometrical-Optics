% EFFICIENCIES_BEAM Scattering of a beam on a spherical particle
%
% Calculation of the trapping efficiencies corresponding to the scattering
% of a focused beam on a spherical particle as a function of the position
% of the particle with respect to the focal point.
%
% See also Ray, BeamGauss, ParticleSpherical.
%
% This example is part of the OTGO - Optical Tweezers in Geometrical Optics
% software package, which complements the article by
% Agnese Callegari, Mite Mijalkov, Burak Gokoz & Giovanni Volpe
% 'Computational toolbox for optical tweezers in geometrical optics'
% (2014).

% Author: Agnese Callegari
% Date: 2014/01/01
% Version: 1.0.0

%% Workspace initialization
clear all;
close all;
clc;

%% Parameters

% Medium
nm = 1.33; % Medium refractive index

% Spherical particle
 Rp = 1e-6; % Particle radius [m] % Particle radius [m]
%  Np = linspace(1,4,8); % Particle refractive index
 Np = 1.50;
% Focusing
f = 100e-6; % Focal length [m]
NA = 1.30; % numerical aperture
L = f*NA/nm; % Iris aperture [m]

% Trapping beam
Ex0 = 1e+4; % x electric field [V/m]
Ey0 = 1i*1e+4; % y electric field [V/m]
w0 = 100e-6; % Beam waist [m]
Nphi = 40; % Azimuthal divisions
Nr = 40; % Radial divisions
power = 5e-3; % Power [W]

dirout= ['.' filesep 'Sphere 40x40 grid f=100' filesep]
mkdir(dirout)
% standard deviation of the spatial variables

% xl = 3e-6; % maximun range along x
% yl = 3e-6; % maximun range along y
% zl = 5e-6; % maximun range along z


% xl = 2e-6; % maximun range along x
% yl = 2e-6; % maximun range along y
% zl = 4e-6; % maximun range along z


%% Initialization

% Trapping beam
bg = BeamGauss(Ex0,Ey0,w0,L,Nphi,Nr);
bg = bg.normalize(power); % Set the power

% Calculates set of rays corresponding to optical beam
r = Ray.beam2focused(bg,f);



% output format : files

% Ndatafiles = length(Rp); % number of files to be generated
%Nsub=length(Np);      % number of rows in each file (elements in a row: x0, y0, z0, fx0, fy0, fz0)
Ntot=102; %  total number of points in the 3D space





Ri=zeros(2,Ntot);
Rpi=zeros(2,Ntot);
xv=zeros(2,Ntot);
yv1=zeros(1,Ntot);



   
    np=Np;
    R=Rp;
    xl = R;
    x=linspace(-4*xl,4*xl,Ntot);
     xs = linspace(-4*xl,0,Ntot/2);
     x = [xs, flip(-xs(1:end-1))];
     y = x;
    bead = ParticleSpherical(Point(0,0,0),R,nm,np);
    Ntot=Ntot-1;
    for i=1:Ntot
         %% X calculation
          %     %Spherical particle
                bead.sp.c.X = x(i);
                bead.sp.c.Y = 0;
                bead.sp.c.Z = 0;
                
                % Calculate force
                forces = bead.force(r);
                force = Vector(x(i),0,0), ...
                    sum(forces.Vx(isfinite(forces.Vx))), ...
                    sum(forces.Vy(isfinite(forces.Vy))), ...
                    sum(forces.Vz(isfinite(forces.Vz))) ...
                    ;
                
%                  fx(i,l,j) = force.Vx;
%                 fxo(i,l,j) = force.Vx;
                fx(j,l) = force.Vx;
                fy(j,l) = force.Vy;
                fz(j,l) = force.Vz;

                
                

    end
        
       
        
        
 