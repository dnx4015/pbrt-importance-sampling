
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


// materials/measuredAdaptive.cpp*
#include "stdafx.h"
#include "materials/measuredAdaptive.h"
#include "paramset.h"
#include "floatfile.h"

AdaptiveHalfangleBRDF::AdaptiveHalfangleBRDF( float *d, 
					uint32_t nth, uint32_t ntd, uint32_t npd, 
					uint32_t mto, uint32_t mpo, uint32_t mth, 
					uint32_t mph, 
					vector<Distribution2DAdaptive*> dist)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), 
		brdf(d),
		nThetaH(nth), nThetaD(ntd), nPhiD(npd), 
		mThetaO(mto), mPhiO(mpo), mThetaH(mth), mPhiH(mph), 
		distribution(dist){
}

Spectrum AdaptiveHalfangleBRDF::f(const Vector &WO, const Vector &WI) 
							const {
	Vector wo = WO, wi = WI;
	Vector wh = wo + wi;
    if (wh.z < 0.f) {
        wo = -wo;
        wi = -wi;
        wh = -wh;
    }
    if (wh.x == 0.f && wh.y == 0.f && wh.z == 0.f) 
		return Spectrum (0.f);
    wh = Normalize(wh);
	
	float whTheta = SphericalTheta(wh);
    float whCosPhi = CosPhi(wh), whSinPhi = SinPhi(wh);
    float whCosTheta = CosTheta(wh), whSinTheta = SinTheta(wh);
    Vector whx(whCosPhi * whCosTheta, whSinPhi * whCosTheta, -whSinTheta);
    Vector why(-whSinPhi, whCosPhi, 0);
    Vector wd(Dot(wi, whx), Dot(wi, why), Dot(wi, wh));

    float wdTheta = SphericalTheta(wd), wdPhi = SphericalPhi(wd);
    if (wdPhi > M_PI) wdPhi -= M_PI;

#define REMAP(V, MAX, COUNT) \
        Clamp(int((V) / (MAX) * (COUNT)), 0, (COUNT)-1)
    int whThetaIndex = REMAP(
						sqrtf(max(0.f, whTheta / (M_PI / 2.f))),
                        1.f, nThetaH);
    int wdThetaIndex = REMAP(wdTheta, M_PI / 2.f, nThetaD);
    int wdPhiIndex = REMAP(wdPhi, M_PI, nPhiD);
#undef REMAP

	int index = wdPhiIndex + nPhiD * 
				(wdThetaIndex + whThetaIndex * nThetaD);
	return Spectrum::FromRGB(&brdf[3*index]);
}

Spectrum AdaptiveHalfangleBRDF::Sample_f(const Vector &wo, 
										Vector *wi, float u1, 
										float u2, float *pdf)const{
	float uv[2], mapPdf;
    float woTheta = SphericalTheta(wo), woPhi = SphericalPhi(wo);
#define REMAP(V, MAX, COUNT) \
	Clamp(int((V) / (MAX) * (COUNT)), 0, (COUNT)-1)
	int woThetaIndex = REMAP(woTheta, (M_PI / 2.f), mThetaO);
	int woPhiIndex = REMAP(woPhi, (M_PI * 2.f), mPhiO);	
#undef REMAP
	int woIndex = woThetaIndex*mPhiO + woPhiIndex;
	distribution[woIndex]->SampleContinuous(u1, u2, uv, &mapPdf);
	if (mapPdf == 0.f) return 0.f;
	
    float theta = uv[1] * (M_PI/2.f), phi = uv[0] * (2.f * M_PI);
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector wh = Normalize(Vector(sintheta * cosphi, sintheta * sinphi, costheta));
    if (!SameHemisphere(wo, wh)) wh = -wh;
	*wi = -wo + 2.f * Dot(wo, wh) * wh;
    
	*pdf = mapPdf / (4.f * Dot(wo, wh));
    if (!SameHemisphere(wo, *wi)) *pdf = 0.f;

    return f(wo, *wi);
}

float AdaptiveHalfangleBRDF::Pdf(const Vector &wo, 
										const Vector &wi) const{
    if (!SameHemisphere(wo, wi)) return 0.f;
	float woTheta = SphericalTheta(wo), woPhi = SphericalPhi(wo);
#define REMAP(V, MAX, COUNT) \
	Clamp(int((V) / (MAX) * (COUNT)), 0, (COUNT)-1)
	int woThetaIndex = REMAP(woTheta, (M_PI / 2.f), mThetaO);
	int woPhiIndex = REMAP(woPhi, (M_PI * 2.f), mPhiO);	
#undef REMAP
	int woIndex = woThetaIndex*mPhiO + woPhiIndex;

	Vector wh = Normalize(wo+wi);
    float whTheta = SphericalTheta(wh), whPhi = SphericalPhi(wh);
	float v = whTheta / (M_PI / 2.f);
	float u = whPhi / (M_PI * 2.f);
	float pdf = distribution[woIndex]->Pdf(u, v) / 
				(4.f * Dot(wo, wh));
	return pdf;
}

// MeasuredAdaptiveMaterial Method Definitions
static map<string, float *> loadedRegularHalfangleAdaptive;
static map<string, vector <Distribution2DAdaptive *> > loadedDistributionAdaptive;

MeasuredAdaptiveMaterial::MeasuredAdaptiveMaterial(
		const string &fName,
      	Reference<Texture<float> > bump, 
		int tp, int mSize, float mPDist, float mRDist) {
    bumpMap = bump;
    const char *suffix = strrchr(fName.c_str(), '.');
    regularHalfangleData = NULL;
    if (!suffix)
        Error("No suffix in measured BRDF filename \"%s\".  "
              "Can't determine file type (.brdf / .merl)", 
			  fName.c_str());
    else {
        // Load RegularHalfangle BRDF Data
       	nThetaH = 90;
        nThetaD = 90;
        nPhiD = 180;
		mThetaO = 32;
		mPhiO = 16;
		mThetaH = 256;
		mPhiH = 32;
        
        if (loadedRegularHalfangleAdaptive.find(fName) != 
			  loadedRegularHalfangleAdaptive.end()) {
            regularHalfangleData = loadedRegularHalfangleAdaptive[fName];
			distribution = loadedDistributionAdaptive[fName];
            return;
        }
        
        FILE *f = fopen(fName.c_str(), "rb");
        if (!f) {
            Error("Unable to open BRDF data file \"%s\"", 
				  fName.c_str());
            return;
        }
        int dims[3];
        if (fread(dims, sizeof(int), 3, f) != 3) {
            Error("Premature end-of-file in measured BRDF data "
				  "file \"%s\"",
                  fName.c_str());
            fclose(f);
            return;
        }
        uint32_t n = dims[0] * dims[1] * dims[2];
        if (n != nThetaH * nThetaD * nPhiD)  {
            Error("Dimensions don't match\n");
            fclose(f);
            return;
        }
        
        regularHalfangleData = new float[3*n];
        const uint32_t chunkSize = 2*nPhiD;
        double *tmp = ALLOCA(double, chunkSize);
        uint32_t nChunks = n / chunkSize;
        Assert((n % chunkSize) == 0);
        float scales[3] = { 1.f/1500.f, 
							1.15f/1500.f, 
							1.66f/1500.f };
        for (int c = 0; c < 3; ++c) {
            int offset = 0;
            for (uint32_t i = 0; i < nChunks; ++i) {
                if (fread(tmp, sizeof(double), chunkSize, f) != 
					  chunkSize) {
                    Error("Premature end-of-file in measured BRDF "
						  "data file \"%s\"",
                          fName.c_str());
                    delete[] regularHalfangleData;
                    regularHalfangleData = NULL;
                    fclose(f);
                    return;
                }
                for (uint32_t j = 0; j < chunkSize; ++j){
                    regularHalfangleData[3 * offset++ + c] = 
						  max(0., tmp[j] * scales[c]);
				}
            }
        }
        loadedRegularHalfangleAdaptive[fName] = regularHalfangleData;
		uint32_t m = mPhiH * mThetaH * mPhiO * mThetaO;
		float* ohHalfangleData = new float[m];
		mapToViewHalfangle(regularHalfangleData, ohHalfangleData);
        fclose(f);
		uint32_t fSize = mThetaH * mPhiH;
		int margSum = 0;
		int condSum = 0;
		Distribution2DAdaptive *temp;
		for (uint32_t i = 0; i < mThetaO; i++){
			for (uint32_t j = 0; j < mPhiO; j++){
				temp = 	new Distribution2DAdaptive(
						&ohHalfangleData[fSize*(j+i*mPhiO)], 
						mPhiH, mThetaH, tp, mSize, mPDist, mRDist);
				distribution.push_back(temp);
				margSum += temp->getMargCount();
				condSum += temp->getCondCount();
			}	
		}
        loadedDistributionAdaptive[fName] = distribution;
		printf("avg theta: %d \n", 
			int(margSum/float(mThetaO*mPhiO)));
		printf("avg phi: %d \n",
			int(condSum/float(mThetaO*mPhiO)));
    }
}

void MeasuredAdaptiveMaterial::mapToViewHalfangle(float* tmpData, 
										float* finalData){
	int index = 0, ind;
	double val, maxVal = 0.0;
	for (uint32_t i = 0; i < mThetaO; i++){
		double thetaOut = i * 0.5 * M_PI / mThetaO;
		for (uint32_t j = 0; j < mPhiO; j++){
			double phiOut = j * 2 * M_PI / mPhiO;
			for (uint32_t k = 0; k < mThetaH; k++){
				double thetaHalf = k * 0.5 * M_PI / mThetaH;
				for (uint32_t l = 0; l < mPhiH; l++){
					double phiHalf = l * 2 * M_PI / mPhiH;
					ind = lookup_brdf_val(tmpData, thetaOut,
												phiOut, thetaHalf,
												phiHalf);
					val = tmpData[ind] + tmpData[ind+1] + 
										tmpData[ind+2];
					finalData[index] = val;
					if (val > maxVal)
						maxVal = val;
					index++;
				}
			}
		}
	}
	/*FILE *f = fopen("function.txt", "w");
	uint32_t m = mPhiH * mThetaH * mPhiO * mThetaO;
	for (uint32_t i = 0; i < m; i++){
		finalData[i] /= maxVal;
		fprintf(f, "%f\n", finalData[i]);
	}
	fclose(f);*/
}

int MeasuredAdaptiveMaterial::lookup_brdf_val(float *brdf, 
										double thetaOut, 
										double phiOut,
										double thetaHalf,
										double phiHalf){
	double hz = cos(thetaHalf);
	double sinThetaHalf = sin(thetaHalf);
	double hx = sinThetaHalf * cos(phiHalf);
	double hy = sinThetaHalf * sin(phiHalf);
	Vector wh(hx, hy, hz);
	wh = Normalize(wh);
	
	double oz = cos(thetaOut);
	double sinThetaOut = sin(thetaOut);
	double ox = sinThetaOut * cos(phiOut);
	double oy = sinThetaOut * sin(phiOut);
	Vector wo(ox, oy, oz);
	wo = Normalize(wo);

    if (!SameHemisphere(wo, wh)) wh = -wh;
	Vector wi = -wo + 2.f * Dot(wo, wh) * wh;
	wi = Normalize(wi);

	float whTheta = SphericalTheta(wh);
    float whCosPhi = CosPhi(wh), whSinPhi = SinPhi(wh);
    float whCosTheta = CosTheta(wh), whSinTheta = SinTheta(wh);
    Vector whx(whCosPhi * whCosTheta, whSinPhi * whCosTheta, 
			   -whSinTheta);
    Vector why(-whSinPhi, whCosPhi, 0);
    Vector wd(Dot(wi, whx), Dot(wi, why), Dot(wi, wh));

    // Compute _index_ into measured BRDF tables
    float wdTheta = SphericalTheta(wd), wdPhi = SphericalPhi(wd);
    if (wdPhi > M_PI) wdPhi -= M_PI;

    // Compute indices _whThetaIndex_, _wdThetaIndex_, _wdPhiIndex_
#define REMAP(V, MAX, COUNT) \
        Clamp(int((V) / (MAX) * (COUNT)), 0, (COUNT)-1)
    int whThetaIndex = REMAP(sqrtf(max(0.f, whTheta / (M_PI / 2.f))),
                             1.f, nThetaH);
    int wdThetaIndex = REMAP(wdTheta, M_PI / 2.f, nThetaD);
    int wdPhiIndex = REMAP(wdPhi, M_PI, nPhiD);
#undef REMAP
    return wdPhiIndex + nPhiD * 
			(wdThetaIndex + whThetaIndex * nThetaD);
}

BSDF *MeasuredAdaptiveMaterial::GetBSDF(
							const DifferentialGeometry &dgGeom,
                            const DifferentialGeometry &dgShading,
                            MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    if (regularHalfangleData)
        bsdf->Add(BSDF_ALLOC(arena, AdaptiveHalfangleBRDF)
            (regularHalfangleData,  
			nThetaH, nThetaD, nPhiD,
			mThetaO, mPhiO, mThetaH, mPhiH, distribution));
    return bsdf;
}


MeasuredAdaptiveMaterial *CreateMeasuredAdaptiveMaterial(
		const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<float> > bumpMap = 
		  mp.GetFloatTextureOrNull("bumpmap");
	int type = mp.FindInt("compression", 0);	
	int maxSize = mp.FindInt("maxSize", 30);
	float minPDist = mp.FindFloat("minPDist", 0.05);
	float minRDist = mp.FindFloat("minRDist", 0.05);
    return new MeasuredAdaptiveMaterial(mp.FindFilename("filename"),
										bumpMap, type, maxSize, 
										minPDist, minRDist);
}
