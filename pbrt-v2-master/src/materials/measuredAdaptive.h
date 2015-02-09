
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

#ifndef PBRT_MATERIALS_MEASURED_ADAPTIVE
#define PBRT_MATERIALS_MEASURED_ADAPTIVE

// materials/measuredAdaptive.h*
#include "pbrt.h"
#include "material.h"
#include "reflection.h"
#include "adaptive.h"
#include "kdtree.h"

class AdaptiveHalfangleBRDF : public BxDF{
public:
	AdaptiveHalfangleBRDF( float *d, 	
					uint32_t nth, uint32_t ntd, uint32_t npd, 
					uint32_t mto, uint32_t mpo, uint32_t mth, 
					uint32_t mph, 
					vector<Distribution2DAdaptive*> dist);
	~AdaptiveHalfangleBRDF(){
		if (brdf)
			delete[] brdf;
		for (uint32_t i = 0; i < distribution.size(); ++i)
			if (distribution[i])
		        delete distribution[i];
	}
	Spectrum f(const Vector &wo, const Vector &wi) const;
	Spectrum Sample_f(const Vector &wo, Vector *wi, 
					float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wi, const Vector &wo) const;
    float *brdf;
    uint32_t nThetaH, nThetaD, nPhiD;
    uint32_t mThetaO, mPhiO, mThetaH, mPhiH;
	vector<Distribution2DAdaptive *> distribution;
};

// MeasuredAdaptiveMaterial Declarations
class MeasuredAdaptiveMaterial : public Material {
public:
    // MeasuredAdaptiveMaterial Public Methods
    MeasuredAdaptiveMaterial(const string &filename, 
						Reference<Texture<float> > bump, int type,
						int mSize, float mPDist, float mRDist);
	~MeasuredAdaptiveMaterial(){
		if (regularHalfangleData)
			delete[]regularHalfangleData;
		for (uint32_t i = 0; i < distribution.size(); ++i)
			if (distribution[i])
		        delete distribution[i];

	}
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
	void mapToViewHalfangle(float *tmpData, float *finalData);
	int lookup_brdf_val(float *brdf, double thetaOut, double phiOut,
						double thetaHalf, double phiHalf);

private:
    // MeasuredAdaptiveMaterial Private Data
    float *regularHalfangleData;
    uint32_t mThetaO, mPhiO, mThetaH, mPhiH;
    uint32_t nThetaH, nThetaD, nPhiD;
	vector<Distribution2DAdaptive *> distribution;
    Reference<Texture<float> > bumpMap;
};


MeasuredAdaptiveMaterial *CreateMeasuredAdaptiveMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_MEASURED_ADAPTIVE
