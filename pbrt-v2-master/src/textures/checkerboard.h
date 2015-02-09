
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_TEXTURES_CHECKERBOARD_H
#define PBRT_TEXTURES_CHECKERBOARD_H

// textures/checkerboard.h*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
#include "montecarlo.h"
#include "shape.h"
#include "parallel.h"
#include "progressreporter.h"

// CheckerboardTexture Declarations
template <typename T> class Checkerboard2DTexture : public Texture<T> {
public:
    // Checkerboard2DTexture Public Methods
    Checkerboard2DTexture(TextureMapping2D *m, Reference<Texture<T> > c1,
                          Reference<Texture<T> > c2, const string &aa)
        : mapping(m), tex1(c1), tex2(c2) {
        // Select antialiasing method for _Checkerboard2DTexture_
        if (aa == "none")             aaMethod = NONE;
        else if (aa == "closedform")  aaMethod = CLOSEDFORM;
        else {
            Warning("Antialiasing mode \"%s\" not understood by "
                    "Checkerboard2DTexture; using \"closedform\"", aa.c_str());
            aaMethod = CLOSEDFORM;
        }
    }
    ~Checkerboard2DTexture() {
        delete mapping;
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        float s, t, dsdx, dtdx, dsdy, dtdy;
        mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
        if (aaMethod == NONE) {
            // Point sample _Checkerboard2DTexture_
            if ((Floor2Int(s) + Floor2Int(t)) % 2 == 0)
                return tex1->Evaluate(dg);
            return tex2->Evaluate(dg);
        }
        else {
            // Compute closed-form box-filtered _Checkerboard2DTexture_ value

            // Evaluate single check if filter is entirely inside one of them
            float ds = max(fabsf(dsdx), fabsf(dsdy));
            float dt = max(fabsf(dtdx), fabsf(dtdy));
            float s0 = s - ds, s1 = s + ds;
            float t0 = t - dt, t1 = t + dt;
            if (Floor2Int(s0) == Floor2Int(s1) && Floor2Int(t0) == Floor2Int(t1)) {
                // Point sample _Checkerboard2DTexture_
                if ((Floor2Int(s) + Floor2Int(t)) % 2 == 0)
                    return tex1->Evaluate(dg);
                return tex2->Evaluate(dg);
            }

            // Apply box filter to checkerboard region
#define BUMPINT(x) \
                (Floor2Int((x)/2) + \
                 2.f * max((x/2)-Floor2Int(x/2) - .5f, 0.f))
            float sint = (BUMPINT(s1) - BUMPINT(s0)) / (2.f * ds);
            float tint = (BUMPINT(t1) - BUMPINT(t0)) / (2.f * dt);
            float area2 = sint + tint - 2.f * sint * tint;
            if (ds > 1.f || dt > 1.f)
                area2 = .5f;
            return (1.f - area2) * tex1->Evaluate(dg) +
                   area2         * tex2->Evaluate(dg);
        }
    }
private:
    // Checkerboard2DTexture Private Data
    TextureMapping2D *mapping;
    Reference<Texture<T> > tex1, tex2;
    enum { NONE, CLOSEDFORM } aaMethod;
};


template <typename T> class Checkerboard3DTexture : public Texture<T> {
public:
    // Checkerboard3DTexture Public Methods
    Checkerboard3DTexture(TextureMapping3D *m, Reference<Texture<T> > c1,
                          Reference<Texture<T> > c2)
        : mapping(m), tex1(c1), tex2(c2) {
    }
    ~Checkerboard3DTexture() {
        delete mapping;
    }
    T Evaluate(const DifferentialGeometry &dg) const {
        Vector dpdx, dpdy;
        Point p = mapping->Map(dg, &dpdx, &dpdy);
        if ((Floor2Int(p.x) + Floor2Int(p.y) + Floor2Int(p.z)) % 2 == 0)
            return tex1->Evaluate(dg);
        else
            return tex2->Evaluate(dg);
    }
private:
    // Checkerboard3DTexture Private Data
    TextureMapping3D *mapping;
    Reference<Texture<T> > tex1, tex2;
};


Texture<float> *CreateCheckerboardFloatTexture(const Transform &tex2world,
        const TextureParams &tp);
Texture<Spectrum> *CreateCheckerboardSpectrumTexture(const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_CHECKERBOARD_H
