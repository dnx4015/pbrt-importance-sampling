#include <set>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <vector>
#include <stdint.h>

#define IOR(x) freopen(x,"r",stdin);
#define FOR(i,a,b) for(int i=a;i<=b;i++)
#define mp make_pair
#define fst first
#define snd second

#define MULT 100000
#define MIN_DIST 0.00001
#define MAX_SIZE 30

using namespace std;

inline float Clamp(float val, float low, float high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}

inline int Clamp(int val, int low, int high) {
    if (val < low) return low;
    else if (val > high) return high;
    else return val;
}

inline int Float2Int(float val) {
    return (int)val;
}

inline float calcDistance(std::pair<int, float>p1, 
						  std::pair<int, float>p2, 
				          std::pair<int, float>p){
	float dist = fabs((p2.fst - p1.fst)*(p1.snd - p.snd)
					- (p1.fst - p.fst)*(p2.snd - p1.snd));
	return dist / sqrt((p2.fst - p1.fst)*(p2.fst - p1.fst)
			+ (p2.snd - p1.snd)*(p2.snd - p1.snd));
}

struct Distribution1D {
    Distribution1D(const float *f, int n) {
        count = n;
        func = new float[n];
		for (int i = 0; i < n; i++)
			func[i] = f[i];
        cdf = new float[n+1];
        cdf[0] = 0.;
        for (int i = 1; i < count+1; ++i)
            cdf[i] = cdf[i-1] + func[i-1] / n;

        funcInt = cdf[count];
        if (funcInt == 0.f) {
            for (int i = 1; i < n+1; ++i)
                cdf[i] = float(i) / float(n);
        }
        else {
            for (int i = 1; i < n+1; ++i){
                cdf[i] /= funcInt;
			}
        }
    }
    ~Distribution1D() {
        delete[] func;
        delete[] cdf;
    }

float SampleContinuous(float u, float *pdf, int *off = NULL) const {
        float *ptr = std::upper_bound(cdf, cdf+count+1, u);
        int offset = max(0, int(ptr-cdf-1));
        if (off) *off = offset;
        float du = (u - cdf[offset]) / (cdf[offset+1] - cdf[offset]);
        if (pdf) *pdf = func[offset] / funcInt;
        return (offset + du) / count;
    }

private:
    friend struct Distribution2D;
    friend struct Distribution1DAdaptive;
    friend struct Distribution2DAdaptive;
    float *func, *cdf;
    float funcInt;
    int count;
};

struct Distribution2D {
    Distribution2D(const float *data, int nu, int nv);
    ~Distribution2D();
    void SampleContinuous(float u0, float u1, float uv[2],
                          float *pdf) const {
        float pdfs[2];
        int v;
        uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
        uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
        *pdf = pdfs[0] * pdfs[1];
    }
    float Pdf(float u, float v) const {
        int iu = Clamp(Float2Int(u * pConditionalV[0]->count), 0,
                       pConditionalV[0]->count-1);
        int iv = Clamp(Float2Int(v * pMarginal->count), 0,
                       pMarginal->count-1);
        if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.f) 
			return 0.f;
        return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
               (pConditionalV[iv]->funcInt * pMarginal->funcInt);
    }
private:
    friend struct Distribution2DAdaptive;
    vector<Distribution1D *> pConditionalV;
    Distribution1D *pMarginal;
};

Distribution2D::Distribution2D(const float *func, int nu, int nv) {
    pConditionalV.reserve(nv);
    for (int v = 0; v < nv; ++v) {
        pConditionalV.push_back(new Distribution1D(&func[v*nu], nu));
    }
    vector<float> marginalFunc;
    marginalFunc.reserve(nv);
    for (int v = 0; v < nv; ++v)
        marginalFunc.push_back(pConditionalV[v]->funcInt);
    pMarginal = new Distribution1D(&marginalFunc[0], nv);
}

Distribution2D::~Distribution2D() {
    delete pMarginal;
    for (uint32_t i = 0; i < pConditionalV.size(); ++i)
        delete pConditionalV[i];
}

struct CDFAdaptive{
	CDFAdaptive(const float *f, int n){
		cdf.insert(mp(0, f[0]));	
		cdf.insert(mp(n, f[n]));	
		float maxDist = 0, dist = 0;
		int index = 0;
		count = 2;
		set<pair<int, float> >::iterator it, nextIt;
		do{
			it = cdf.begin();
			nextIt = cdf.begin();
			nextIt++;
			maxDist = 0;
			for(nextIt, it; nextIt != cdf.end(); 
				  it++, nextIt++){
				FOR(i, (it->fst)+1, (nextIt->fst)-1){
					dist = calcDistance(
						mp(it->fst, it->snd*MULT), 
						mp(nextIt->fst, nextIt->snd*MULT), 
						mp(i, f[i] * MULT));
					if (dist > maxDist){
						maxDist = dist;
						index = i;
					}
				}	
			}
			if (maxDist > MIN_DIST){
				cdf.insert(mp(index, f[index]));
				count += 1;
			}
		}while(maxDist > MIN_DIST && cdf.size() < MAX_SIZE);
	}
	
	void printCdfAdap(){
		set<pair<int, float> >::iterator it = cdf.begin(), 
		nextIt = cdf.begin();
		nextIt++;
		for(nextIt, it; nextIt != cdf.end(); it++, nextIt++)
			printf("%d, %f -> %d, %f\n",it->fst, it->snd, 
					nextIt->fst, nextIt->snd);
		
	}
	
	void toArrays(float *cdfVal, int *cdfIndex){
		set<pair<int, float> >::iterator it;
		int i = 0;
		for(it = cdf.begin(); it != cdf.end(); it++){
			cdfVal[i] = it->snd;
			cdfIndex[i] = it->fst;
			i++;
		}
	}

	void updateFunc(const float *f, float *func){
		set<pair<int, float> >::iterator it = cdf.begin(), 
		nextIt = cdf.begin();
		nextIt++;
		int index1, index2, i = 0, n;
		float sumFunc;
		for(nextIt, it; nextIt != cdf.end(); it++, nextIt++){
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

private:
	friend struct Distribution1DAdaptive;
	set< pair<int, float> > cdf;
 	int count;
};

// Monte Carlo Utility Declarations
struct Distribution1DAdaptive {
	Distribution1DAdaptive(Distribution1D *pMarginal) {
        uniformCount = pMarginal->count;
		funcInt = pMarginal->funcInt;
       	CDFAdaptive *cdfAdaptive = new CDFAdaptive(pMarginal->cdf, uniformCount);
		adaptiveCount = cdfAdaptive->count;
		func = new float[adaptiveCount];
		cdf = new float[adaptiveCount];
		index = new int[adaptiveCount];
		cdfAdaptive->toArrays(cdf, index);
		cdfAdaptive->updateFunc(pMarginal->func, func);
		//cdfAdaptive->printCdfAdap();
    }

    ~Distribution1DAdaptive() {
        delete[] cdf;
        delete[] func;
        delete[] index;
    }
    
	float SampleContinuous(float u, float *pdf, int *off = NULL) 
		const {
	    float *ptr = std::upper_bound(cdf, cdf+adaptiveCount+1, u);
        int offset = max(0, int(ptr-cdf-1));
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
	
	float Pdf(float u, int &adaptiveU){
		int offset = Clamp(Float2Int(u * uniformCount), 0,uniformCount-1);
		if (funcInt == 0)
			return 0.f;
		int *ptr = std::upper_bound(index, index+adaptiveCount+1, offset);
		adaptiveU = max(0, int(ptr - index - 1));
        return func[adaptiveU] / funcInt;
	}

private:
    friend struct Distribution2DAdaptive;
	//CDFAdaptive *cdfAdaptive;
	int uniformCount;
	int adaptiveCount;
	float *cdf, *func;
	int *index;
	float funcInt;
};

struct Distribution2DAdaptive {
    Distribution2DAdaptive(Distribution2D *dist);
    ~Distribution2DAdaptive();
    void SampleContinuous(float u0, float u1, float uv[2],
                          float *pdf) const;
	float Pdf(float u, float v) const;
private:
    vector<Distribution1DAdaptive *> pConditionalV;
    Distribution1DAdaptive *pMarginal;
};

Distribution2DAdaptive::Distribution2DAdaptive(Distribution2D *dist) {
	pMarginal = new Distribution1DAdaptive(dist->pMarginal);
	int nu = dist->pConditionalV[0]->count;
	int nv = dist->pMarginal->count;
	float *func = new float[nu];
	//int w = 1;
	int na = pMarginal->adaptiveCount - 1;
	for(int k = 0; k < na; k++){
		for (int i = 0; i < nu; i++){
			for (int j = pMarginal->index[k]; 
  					j < pMarginal->index[k+1]; j++){
				func[i] = dist->pConditionalV[j]->func[i];
			}
		}
		//printf("\n\n%d\n\n", w++);
		pConditionalV.push_back(
			new Distribution1DAdaptive(
				new Distribution1D(func, na)));
	}
}

void Distribution2DAdaptive::SampleContinuous(float u0, float u1, 
  											  float uv[2], float *pdf) 
const {
	float pdfs[2];
	int v;
	uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
	uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
	*pdf = pdfs[0] * pdfs[1];
}
    
float Distribution2DAdaptive::Pdf(float u, float v) const {
	int adaptiveV = 0, adaptiveU = 0;
	float pdfV = pMarginal->Pdf(v, adaptiveV);
	float pdfU = pConditionalV[adaptiveV]->Pdf(u, adaptiveU);
	return pdfV * pdfU;
}

int main(){
	IOR("function.in");
	int n = 32;
	int m = 256;
	float *func = new float[n*m];
	for (int i = 0; i < n*m; i++){
		scanf("%f", &func[i]);	
	}
	Distribution2D *dist = new Distribution2D(func, n, m);
	Distribution2DAdaptive *distAdap = 
		new Distribution2DAdaptive(dist);
	float uv[2];
	float pdf;
	distAdap->SampleContinuous(0.4, 0.7, uv, &pdf);
	printf("u: %f, v: %f\n", uv[0], uv[1]);
	printf("pdf:%f \n", pdf);
}
