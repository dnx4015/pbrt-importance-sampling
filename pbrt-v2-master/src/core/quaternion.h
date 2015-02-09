
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

#ifndef PBRT_CORE_QUATERNION_H
#define PBRT_CORE_QUATERNION_H

// core/quaternion.h*
#include "pbrt.h"
#include "geometry.h"

// Quaternion Declarations
struct Quaternion {
    // Quaternion Public Methods
    Quaternion() { v = Vector(0., 0., 0.); w = 1.f; }
    Quaternion &operator+=(const Quaternion &q) {
        v += q.v;
        w += q.w;
        return *this;
    }
    friend Quaternion operator+(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret += q2;
    }
    Quaternion &operator-=(const Quaternion &q) {
        v -= q.v;
        w -= q.w;
        return *this;
    }
    friend Quaternion operator-(const Quaternion &q1, const Quaternion &q2) {
        Quaternion ret = q1;
        return ret -= q2;
    }
    Quaternion &operator*=(float f) {
        v *= f;
        w *= f;
        return *this;
    }
    Quaternion operator*(float f) const {
        Quaternion ret = *this;
        ret.v *= f;
        ret.w *= f;
        return ret;
    }
    Quaternion &operator/=(float f) {
        v /= f;
        w /= f;
        return *this;
    }
    Quaternion operator/(float f) const {
        Quaternion ret = *this;
        ret.v /= f;
        ret.w /= f;
        return ret;
    }
    Transform ToTransform() const;
    Quaternion(const Transform &t);

    // Quaternion Public Data
    Vector v;
    float w;
};


Quaternion Slerp(float t, const Quaternion &q1, const Quaternion &q2);

// Quaternion Inline Functions
inline Quaternion operator*(float f, const Quaternion &q) {
    return q * f;
}


inline float Dot(const Quaternion &q1, const Quaternion &q2) {
    return Dot(q1.v, q2.v) + q1.w * q2.w;
}


inline Quaternion Normalize(const Quaternion &q) {
    return q / sqrtf(Dot(q, q));
}



#endif // PBRT_CORE_QUATERNION_H
