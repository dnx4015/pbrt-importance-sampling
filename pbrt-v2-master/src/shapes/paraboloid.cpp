
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/paraboloid.cpp*
#include "stdafx.h"
#include "shapes/paraboloid.h"
#include "paramset.h"

// Paraboloid Method Definitions
Paraboloid::Paraboloid(const Transform *o2w, const Transform *w2o, bool ro,
                       float rad, float z0, float z1,
                       float tm)
    : Shape(o2w, w2o, ro) {
    radius = rad;
    zmin = min(z0,z1);
    zmax = max(z0,z1);
    phiMax = Radians(Clamp(tm, 0.0f, 360.0f));
}


BBox Paraboloid::ObjectBound() const {
    Point p1 = Point( -radius, -radius, zmin );
    Point p2 = Point(  radius,  radius, zmax );
    return BBox( p1, p2 );
}


bool Paraboloid::Intersect(const Ray &r, float *tHit,
        float *rayEpsilon, DifferentialGeometry *dg) const {
    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute quadratic paraboloid coefficients
    float k = zmax/(radius*radius);
    float A =   k*(ray.d.x * ray.d.x + ray.d.y * ray.d.y);
    float B = 2*k*(ray.d.x * ray.o.x + ray.d.y * ray.o.y) -
                ray.d.z;
    float C =   k*(ray.o.x * ray.o.x + ray.o.y * ray.o.y) -
                ray.o.z;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute paraboloid inverse mapping
    phit = ray(thit);
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;

    // Test paraboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        if (thit == t1) return false;
        thit = t1;
        if (t1 > ray.maxt) return false;
        // Compute paraboloid inverse mapping
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        if (phit.z < zmin || phit.z > zmax || phi > phiMax)
            return false;
    }

    // Find parametric representation of paraboloid hit
    float u = phi / phiMax;
    float v = (phit.z-zmin) / (zmax-zmin);

    // Compute parabaloid $\dpdu$ and $\dpdv$
    Vector dpdu(-phiMax * phit.y, phiMax * phit.x, 0.);
    Vector dpdv = (zmax - zmin) *
        Vector(phit.x / (2.f * phit.z), phit.y / (2.f * phit.z), 1.);

    // Compute parabaloid $\dndu$ and $\dndv$
    Vector d2Pduu = -phiMax * phiMax *
                    Vector(phit.x, phit.y, 0);
    Vector d2Pduv = (zmax - zmin) * phiMax *
                    Vector(-phit.y / (2.f * phit.z),
                           phit.x / (2.f * phit.z),
                           0);
    Vector d2Pdvv = -(zmax - zmin) * (zmax - zmin) *
                    Vector(phit.x / (4.f * phit.z * phit.z),
                           phit.y / (4.f * phit.z * phit.z),
                           0.);

    // Compute coefficients for fundamental forms
    float E = Dot(dpdu, dpdu);
    float F = Dot(dpdu, dpdv);
    float G = Dot(dpdv, dpdv);
    Vector N = Normalize(Cross(dpdu, dpdv));
    float e = Dot(N, d2Pduu);
    float f = Dot(N, d2Pduv);
    float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    float invEGF2 = 1.f / (E*G - F*F);
    Normal dndu = Normal((f*F - e*G) * invEGF2 * dpdu +
                         (e*F - f*E) * invEGF2 * dpdv);
    Normal dndv = Normal((g*F - f*G) * invEGF2 * dpdu +
                         (f*F - g*E) * invEGF2 * dpdv);

    // Initialize _DifferentialGeometry_ from parametric information
    const Transform &o2w = *ObjectToWorld;
    *dg = DifferentialGeometry(o2w(phit), o2w(dpdu), o2w(dpdv),
                               o2w(dndu), o2w(dndv), u, v, this);

    // Update _tHit_ for quadric intersection
    *tHit = thit;

    // Compute _rayEpsilon_ for quadric intersection
    *rayEpsilon = 5e-4f * *tHit;
    return true;
}



bool Paraboloid::IntersectP(const Ray &r) const {
    float phi;
    Point phit;
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

    // Compute quadratic paraboloid coefficients
    float k = zmax/(radius*radius);
    float A =   k*(ray.d.x * ray.d.x + ray.d.y * ray.d.y);
    float B = 2*k*(ray.d.x * ray.o.x + ray.d.y * ray.o.y) -
                ray.d.z;
    float C =   k*(ray.o.x * ray.o.x + ray.o.y * ray.o.y) -
                ray.o.z;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    // Compute intersection distance along ray
    if (t0 > ray.maxt || t1 < ray.mint)
        return false;
    float thit = t0;
    if (t0 < ray.mint) {
        thit = t1;
        if (thit > ray.maxt) return false;
    }

    // Compute paraboloid inverse mapping
    phit = ray(thit);
    phi = atan2f(phit.y, phit.x);
    if (phi < 0.) phi += 2.f*M_PI;

    // Test paraboloid intersection against clipping parameters
    if (phit.z < zmin || phit.z > zmax || phi > phiMax) {
        if (thit == t1) return false;
        thit = t1;
        if (t1 > ray.maxt) return false;
        // Compute paraboloid inverse mapping
        phit = ray(thit);
        phi = atan2f(phit.y, phit.x);
        if (phi < 0.) phi += 2.f*M_PI;
        if (phit.z < zmin || phit.z > zmax || phi > phiMax)
            return false;
    }
    return true;
}


float Paraboloid::Area() const {
    return phiMax/12.0f *
        (powf(1+4*zmin, 1.5f) - powf(1+4*zmax, 1.5f));
}


Paraboloid *CreateParaboloidShape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    float radius = params.FindOneFloat("radius", 1);
    float zmin = params.FindOneFloat("zmin", 0);
    float zmax = params.FindOneFloat("zmax", 1);
    float phimax = params.FindOneFloat("phimax", 360);
    return new Paraboloid(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax);
}


