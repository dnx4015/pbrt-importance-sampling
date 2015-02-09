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


// core/adaptive.cpp*
#include "stdafx.h"
#include "adaptive.h"
#include "montecarlo.h"

enum Compression {DP, VW, RD, PD, RW, OP, LA};

void CDFAdaptive::DouglasPeucker(const float *f, int n){
	cdf.insert(mp(0, f[0]));
	cdf.insert(mp(n, f[n]));
	float maxDist = 0, dist = 0;
	count = 2;
	set<pair<int, float> >::iterator it, nextIt;
	int ind = 0;
	do{
		it = cdf.begin();
		nextIt = cdf.begin();
		nextIt++;
		maxDist = -100;
		for(; nextIt != cdf.end(); 
			  it++, nextIt++){
			FOR(i, (it->fst)+1, (nextIt->fst)-1){
				dist = calcDistance(
					mp(it->fst, it->snd*MULT), 
					mp(nextIt->fst, nextIt->snd*MULT), 
					mp(i, f[i] * MULT));
				if (dist > maxDist){
					maxDist = dist;
					ind = i;
				}
			}	
		}
		if (maxDist > minPDist){
			cdf.insert(mp(ind, f[ind]));
			count += 1;
		}
	}while(maxDist > minPDist && int(cdf.size()) < maxSize);
}

int CDFAdaptive::trianglArea(int i, int j, int k, const float *f){
	Point a = Point(float(i), f[i], 1.0);
	Point b = Point(float(j), f[j], 1.0);
	Point c = Point(float(k), f[k], 1.0);
	Vector ab = b - a;
	Vector ac = c - a;
	Vector result = Cross(ab, ac);
	return int(result.Length() * 0.5 * MULT);
}

void CDFAdaptive::VisvalingamWhyatt(const float *f, int n){
	int numRemove = n - maxSize, minIndex, minArea, area;
	vector<int> temp;
	FOR(i, 0, n)
		temp.pb(i);
	FOR(i, 0, numRemove){
		minIndex = 1;
		minArea = trianglArea(temp[0], temp[1], temp[2], f);
		FOR(j, 2, int(temp.size() - 2)){
			area = trianglArea(temp[j-1], temp[j], temp[j+1], f);
			if (area < minArea){
				minArea = area;
				minIndex = j;
			}
		}
		temp.erase(temp.begin() + minIndex);
	}
	
	sort(temp.begin(), temp.end());

	FOR(i, 0, maxSize-1){
		cdf.insert(mp(temp[i], f[temp[i]]));
		count++;
	}
}

float distance(int i, int j, const float *f){
	Point a = Point(float(i), f[i], 1.0);
	Point b = Point(float(j), f[j], 1.0);
	Vector dist = a-b;
	return floor(dist.Length() * MULT);	
}

void CDFAdaptive::RadialDistance(const float *f, int n){
	int key = 0, test = 1;
	cdf.insert(mp(key, f[key]));
	count = 1;
	float d;
	while(key <= n){
		test = key+1;
		d = distance(key, test, f);
		while (test < n && d < minRDist){
			test++;
			d = distance(key, test, f);
		}
		key = test;
		if(key <= n){
			cdf.insert(mp(key, f[key]));
			count++;
		}
	}
}

void CDFAdaptive::PerpendicularDistance(const float *f, int n){
	vector<int> temp;
	FOR(i, 0, n)
		temp.pb(i);
	float dist;
	FOR(i, 1, int(temp.size())-2){
		dist = calcDistance(
				mp(i-1, f[i-1]*MULT), 
				mp(i+1, f[i+1]*MULT), 
				mp(i, f[i] * MULT));
		if (dist < minPDist){
			temp.erase(temp.begin() + i);
			i--;
		}
	}
	FOR(i, 0, int(temp.size() - 1)){
		cdf.insert(mp(temp[i], f[temp[i]]));
		count++;
	}
}


void CDFAdaptive::ReumannWitkam(const float *f, int n){
	vector<int> temp;
	FOR(i, 0, n)
		temp.pb(i);
	int key = 0, aux = 1, test = 2;
	float dist;
	while(key <= int(temp.size()) - 3){
		aux = key + 1;
		test = key + 2;
		dist = calcDistance(mp(key, f[key]*MULT), 
							mp(aux, f[aux]*MULT), 
							mp(test, f[test] * MULT));
		while (test < n && dist < minPDist){
			test++;
			dist = calcDistance(mp(key, f[key]*MULT), 
								mp(aux, f[aux]*MULT), 
								mp(test, f[test] * MULT));
		}
		if(key + 1 > test - 1)
			temp.erase(temp.begin()+key+1, temp.begin()+test-1);
		key++;
	}

	FOR(i, 0, int(temp.size())-1){
		cdf.insert(mp(temp[i], f[temp[i]]));
		count++;
	}
}

void CDFAdaptive::Opheim(const float *f, int n){
	vector<int> temp;
	FOR(i, 0, n)
		temp.pb(i);
	int key = 0, aux = 1, test = 2;
	float pDist, rDist;
	while(key <= int(temp.size()) - 3){
		aux = key + 1;
		test = key + 2;
		pDist = calcDistance(mp(key, f[key]*MULT), 
							mp(aux, f[aux]*MULT), 
							mp(test, f[test] * MULT));
		while (test < n && pDist < minPDist && rDist < minRDist){
			test++;
			pDist = calcDistance(mp(key, f[key]*MULT), 
								mp(aux, f[aux]*MULT), 
								mp(test, f[test] * MULT));
			rDist = distance(key, test, f);
		}
		temp.erase(temp.begin()+key+1, temp.begin()+test-1);
		key++;
	}

	FOR(i, 0, int(temp.size())-1){
		cdf.insert(mp(temp[i], f[temp[i]]));
		count++;
	}
}

void CDFAdaptive::Lang(const float *f, int n){
	vector<int> temp;
	FOR(i, 0, n)
		temp.pb(i);
	int key = 0, end = n;
	int maxIndex = key+1;
	float dist, maxDist;
	do{
		if (key + 1 == end){
			if (end != int(temp.size()-1)){
				key = end;
				end = int(temp.size()-1);
			}
		}else{
			maxIndex = key + 1;
			dist = calcDistance(mp(key, f[key]*MULT), 
								mp(end, f[end]*MULT), 
								mp(maxIndex, f[maxIndex] * MULT));
			maxDist = dist;
			FOR(i, maxIndex+1, end-1){
				dist = calcDistance(mp(key, f[key]*MULT), 
									mp(end, f[end]*MULT), 
									mp(i, f[i] * MULT));
				if (dist > maxDist){
					maxDist = dist;
					maxIndex =  i;
				}
			}
			
			if(maxDist > minPDist){
				end--;	
			}else{
				temp.erase(temp.begin()+key+1, temp.begin()+end);
				key = end;
				end = int(temp.size()-1);
			}
		}
	} while (key < int(temp.size()-2) || end != int(temp.size()-1));

	FOR(i, 0, int(temp.size())-1){
		cdf.insert(mp(temp[i], f[temp[i]]));
		count++;
	}
}

CDFAdaptive::CDFAdaptive(const float *f, int n, int tp, int mSize, 
						float mDist, float mRDist){
	Compression type = Compression(tp);
	maxSize = mSize;
	minPDist = mDist;
	minRDist = mRDist;
	count = 0;
	switch (type){
		case DP:
			DouglasPeucker(f, n);
			break;
		case VW:
			VisvalingamWhyatt(f, n);
			break;
		case RD:
			RadialDistance(f, n);
			break;
		case PD:
			PerpendicularDistance(f, n);
			break;
		case RW:
			ReumannWitkam(f, n);
			break;
		case OP:
			Opheim(f, n);
			break;
		case LA:
			Lang(f, n);
			break;
	}
}

void CDFAdaptive::printfCdfAdap(){
	set<pair<int, float> >::iterator it = cdf.begin(), 
	nextIt = cdf.begin();
	nextIt++;
	for(; nextIt != cdf.end(); it++, nextIt++)
		printf("%d, %f -> %d, %f\n",it->fst, it->snd, 
				nextIt->fst, nextIt->snd);
	
}

void CDFAdaptive::toArrays(float *cdfVal, int *cdfIndex){
	set<pair<int, float> >::iterator it;
	int i = 0;
	for(it = cdf.begin(); it != cdf.end(); it++){
		cdfVal[i] = it->snd;
		cdfIndex[i] = it->fst;
		i++;
	}
}

void CDFAdaptive::updateFunc(const float *f, float *func){
	set<pair<int, float> >::iterator it = cdf.begin(), 
	nextIt = cdf.begin();
	nextIt++;
	int index1, index2, i = 0, n;
	float sumFunc;
	for(; nextIt != cdf.end(); it++, nextIt++){
		index1 = it->fst;
		index2 = nextIt->fst;
		n = index2-index1;
		sumFunc = 0.0;
		FOR(j, index1, index2-1){
			sumFunc += f[j];
		}
		func[i++] = sumFunc/n;
	}
}

Distribution1DAdaptive::Distribution1DAdaptive(const float *f, 
		int n, int type, int maxSize, float minPDist, 
		float minRDist){
	Distribution1D *dist = new Distribution1D(f, n);
	uniformCount = dist->count;
	funcInt = dist->funcInt;
	CDFAdaptive *cdfAdaptive = new CDFAdaptive(dist->cdf, 
											   uniformCount, type,
											   maxSize, 
											   minPDist, minRDist);
	adaptiveCount = cdfAdaptive->count;
	func = new float[adaptiveCount];
	cdf = new float[adaptiveCount];
	index = new int[adaptiveCount];
	cdfAdaptive->toArrays(cdf, index);
	cdfAdaptive->updateFunc(dist->func, func);
}

Distribution1DAdaptive::Distribution1DAdaptive(
		Distribution1D *dist, int type, int maxSize, 
		float minPDist,	float minRDist) {
	uniformCount = dist->count;
	funcInt = dist->funcInt;
	CDFAdaptive *cdfAdaptive = new CDFAdaptive(dist->cdf, 
											   uniformCount, type, 
											   maxSize, minPDist,
											   minRDist);
	adaptiveCount = cdfAdaptive->count;
	func = new float[adaptiveCount];
	cdf = new float[adaptiveCount];
	index = new int[adaptiveCount];
	cdfAdaptive->toArrays(cdf, index);
	cdfAdaptive->updateFunc(dist->func, func);
}

float Distribution1DAdaptive::SampleContinuous(float u, float *pdf, 
											   int *off = NULL) {
	float *ptr = std::upper_bound(cdf, cdf+adaptiveCount, u);
	int offset = max(0, int(ptr-cdf-1));
	offset = min(adaptiveCount-1, offset);
	if (off) *off = offset;
	float du = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
	if (pdf){
		if (funcInt != 0.f)
			*pdf = func[offset] / funcInt;
		else
			*pdf = 0.f;
	}
	return (index[offset] + du) / uniformCount;
}

float Distribution1DAdaptive::Pdf(float u, int &adaptiveU){
	int offset = Clamp(Float2Int(u * uniformCount), 
						0, uniformCount-1);
	if (funcInt == 0)
		return 0.f;
	int *ptr = std::upper_bound(index, index+adaptiveCount, 
								offset);
	adaptiveU = max(0, int(ptr - index - 1));
	adaptiveU = min(adaptiveCount-1, adaptiveU);
	return func[adaptiveU] / funcInt;
}

Distribution2DAdaptive::Distribution2DAdaptive(const float *data, 			int nu, int nv, int type, int maxSize, 
		float minPDist, float minRDist) {
	Distribution2D *dist = new Distribution2D(data, nu, nv);
	pMarginal = new Distribution1DAdaptive(dist->pMarginal, type,
											maxSize, minPDist,
											minRDist);
	Distribution1DAdaptive *temp;
	float *func = new float[nu];
	int na = pMarginal->adaptiveCount - 1;
	int sum = 0;
	for(int k = 0; k < na; k++){
		for (int i = 0; i < nu; i++){
			func[i] = 0;
			for (int j = pMarginal->index[k]; 
  					j < pMarginal->index[k+1]; j++){
				func[i] += dist->pConditionalV[j]->func[i];
			}
			func[i] /= pMarginal->index[k+1] - pMarginal->index[k];
		}
		temp = new Distribution1DAdaptive(func, nu, type, maxSize, 
										minPDist, minRDist);
		pConditionalV.push_back(temp);
		sum += (temp->adaptiveCount - 1);
	}
	marginalCount = na;
	conditionalCount = int(sum / float(na));
}

void Distribution2DAdaptive::SampleContinuous(float u0, float u1, 
  											  float uv[2], 
											  float *pdf) {
	float pdfs[2];
	int v;
	uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
	uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
	*pdf = pdfs[0] * pdfs[1];
}
    
float Distribution2DAdaptive::Pdf(float u, float v) {
	int adaptiveV = 0, adaptiveU = 0;
	float pdfV = pMarginal->Pdf(v, adaptiveV);
	float pdfU = pConditionalV[adaptiveV]->Pdf(u, adaptiveU);
	return pdfV * pdfU;
}

Distribution2DAdaptive::~Distribution2DAdaptive() {
//	if (pMarginal)
	    delete pMarginal;
    for (uint32_t i = 0; i < pConditionalV.size(); ++i)
//		if(pConditionalV[i])
	        delete pConditionalV[i];
}


