#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>

#define REPEAT(token, num) for (token = 0; token < num; token++)
static inline int64_t min(int64_t length1, int64_t length2) { return (length1 < length2) ? length1 : length2; }
static inline int64_t max(int64_t length1, int64_t length2) { return (length1 > length2) ? length1 : length2; }

//#define DEBUG
#ifdef DEBUG
#define PRINTF printf
#else
#define PRINTF(...)
#endif

typedef int64_t num;
typedef int64_t num_nums;
typedef int64_t test_cases;

num calcModuloExp(num base, num exp, num MOD) {
    num result = 1, cur = base;
    while (exp) {
        if (exp & 1) result *= cur, result %= MOD;
        cur *= cur, cur %= MOD, exp >>= 1;
    }
    return result;
}

num calcModuloInverse(num dividend, num MOD) {
    num lastRemainder = dividend, lastInverse = 1, remainder = MOD, inverse = 0, quotient, temp;
    while (remainder) {
        quotient = lastRemainder/remainder;
        temp = remainder, remainder = lastRemainder % remainder, lastRemainder = temp;
        temp = inverse, inverse = lastInverse-quotient*inverse, lastInverse = temp;
    }
    if (lastInverse < 0) lastInverse += MOD;
    return lastInverse;
}

num MOD1 = 998244353, rootsOfUnity1[100000], r1 = 119, k1 = 23, g1 = 3;
num MOD2 = 2013265921, rootsOfUnity2[100000], r2 = 15, k2 = 27, g2 = 31;
void calcRootsOfUnity(num MOD, num *rootsOfUnity, num r, num k, num g, num arrLengthPower) {
    num rootOfUnity = calcModuloExp(g, r << (k-arrLengthPower), MOD);
    num_nums i;
    rootsOfUnity[0] = 1;
    REPEAT(i, 1 << arrLengthPower) rootsOfUnity[i+1] = (rootsOfUnity[i]*rootOfUnity) % MOD;
}

num bitReversed[100000];
void computeBitReversals(num power) {
    num_nums max = 1 << power, i;
    REPEAT(i, max) bitReversed[i] = (bitReversed[i >> 1] | ((i & 1) << power)) >> 1;
}

void NTT(num MOD, num *rootsOfUnity, num_nums maxArrLengthPower, num *arr, num_nums arrLengthPower, num *ntt, bool inverse) {
    if (arrLengthPower == 0) {
        ntt[0] = arr[0];
        return;
    }
    num_nums i, arrLength = 1 << arrLengthPower, maxArrLength = 1 << maxArrLengthPower;
    num storeNum;
    memcpy(ntt, arr, arrLength*sizeof(num));
    REPEAT(i, arrLength) if (i < bitReversed[i]) {
        storeNum = ntt[i];
        ntt[i] = ntt[bitReversed[i]];
        ntt[bitReversed[i]] = storeNum;
    }
    
    num_nums j, kPower, k, halfK, jumpSizePower;
    //In each iteration, the FFTs of groups of k elements are computed.
    //(i.e. groups of 2, then groups of 4, then groups of 8, etc.)
    for (int kPower = 1; kPower <= arrLengthPower; kPower++) {
        k = 1 << kPower;
        halfK = 1 << (kPower-1);
        jumpSizePower = maxArrLengthPower-kPower;
        //In each iteration, FFT of one group of k elements is computed.
        for (int j = 0; j < arrLength; j += k) {
            //In each iteration, FFT values for ntt[i+j] and ntt[i+j+halfK] are computed.
            for (int i = 0; i < halfK; i++) {
                //Let ind be either i+j or i+j+halfK
                //FFT_k[ind] = FFT_even,halfK[ind % halfK] + rootOfUnity*FFT_odd,halfK[ind % halfK]
                //FFT_even,halfK[ind] is old value of ntt[i+j]
                //FFT_odd,halfK[ind % halfK] is old value of ntt[i+j+halfK]
                //rootOfUnity is positive for ind=i+j, negative for ind=i+j+halfK
                storeNum = (rootsOfUnity[inverse ? (maxArrLength-(i << jumpSizePower)) : (i << jumpSizePower)]*ntt[i+j+halfK]) % MOD;
                //Note that order of operations is crucial: Otherwise, ntt[i+j] will be modified before it can be used to compute ntt[i+j+halfK]
                ntt[i+j+halfK] = ntt[i+j]-storeNum;
                if (ntt[i+j+halfK] < 0) ntt[i+j+halfK] += MOD;
                ntt[i+j] += storeNum;
                if (ntt[i+j] > MOD) ntt[i+j] -= MOD;
            }
        }
    }
    if (inverse) {
        num_nums inverseLength = calcModuloInverse(arrLength, MOD);
        REPEAT(i, arrLength) ntt[i] *= inverseLength, ntt[i] %= MOD;
    }
}

enum { ARRAY_SIZE = 32768, ARRAY_SIZE_POWER = 15 };
num poly[ARRAY_SIZE], poly2[ARRAY_SIZE], poly3a[ARRAY_SIZE], poly3b[ARRAY_SIZE], poly3[ARRAY_SIZE], ntt1a[ARRAY_SIZE], ntt2a[ARRAY_SIZE], ntt3a[ARRAY_SIZE], ntt1b[ARRAY_SIZE], ntt2b[ARRAY_SIZE], ntt3b[ARRAY_SIZE];
num_nums degree;

int main() {
    freopen("input.txt","r",stdin);
    num_nums i;
    test_cases numTestCases = 1, l;
    num multiplier, mod1Inverse;
    mod1Inverse = calcModuloInverse(MOD1, MOD2);
    calcRootsOfUnity(MOD1, rootsOfUnity1, r1, k1, g1, ARRAY_SIZE_POWER);
    calcRootsOfUnity(MOD2, rootsOfUnity2, r2, k2, g2, ARRAY_SIZE_POWER);
    computeBitReversals(ARRAY_SIZE_POWER);
    
    REPEAT(l, numTestCases) {
        scanf("%" PRId64, &degree);
        REPEAT(i, degree+1) scanf("%" PRId64, poly+i);
        REPEAT(i, degree+1) scanf("%" PRId64, poly2+i);
        for (i = degree+1; i < ARRAY_SIZE; i++) poly[i] = poly2[i] = 0;
        
        NTT(MOD1, rootsOfUnity1, ARRAY_SIZE_POWER, poly, ARRAY_SIZE_POWER, ntt1a, false);
        NTT(MOD1, rootsOfUnity1, ARRAY_SIZE_POWER, poly2, ARRAY_SIZE_POWER, ntt2a, false);
        REPEAT(i, ARRAY_SIZE) ntt3a[i] = (ntt1a[i]*ntt2a[i]) % MOD1;
        NTT(MOD1, rootsOfUnity1, ARRAY_SIZE_POWER, ntt3a, ARRAY_SIZE_POWER, poly3a, true);
        
        NTT(MOD2, rootsOfUnity2, ARRAY_SIZE_POWER, poly, ARRAY_SIZE_POWER, ntt1b, false);
        NTT(MOD2, rootsOfUnity2, ARRAY_SIZE_POWER, poly2, ARRAY_SIZE_POWER, ntt2b, false);
        REPEAT(i, ARRAY_SIZE) ntt3b[i] = (ntt1b[i]*ntt2b[i]) % MOD2;
        NTT(MOD2, rootsOfUnity2, ARRAY_SIZE_POWER, ntt3b, ARRAY_SIZE_POWER, poly3b, true);
        
        REPEAT(i, ARRAY_SIZE) {
            multiplier = mod1Inverse*(poly3b[i]-poly3a[i]);
            multiplier %= MOD2;
            if (multiplier < 0) multiplier += MOD2;
            poly3[i] = poly3a[i]+multiplier*MOD1;
        }
        
        REPEAT(i, 2*degree+1) printf("%" PRId64 "%c", poly3[i], (i == 2*degree) ? '\n' : ' ');
    }
    
    exit(0);
}