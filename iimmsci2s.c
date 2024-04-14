/* iimmsci2s.c
 * 
 * Returns parameters corresponding to the MSCI given the MSCM model, 
 * or vice versa, for the case of 2 species by minimising the KL divergence.
 *
 * Adapted from IMMSci2s.c by Ziheng Yang; see Jiao et al (2020, Syst. Biol.)
 * 
 * Yuttapong Thawornwattana (2024)
 *
 * cl -Ox iimmsci2s.c tools.c
 * icc -o iimmsci2s -O3 iimmsci2s.c tools.c -lm
 *
 * Usage: iimmsci2s mode out klt_opt im_model n
 *
 * Model descriptions:
 * ## MSCI model ##
 * ((A)S:tauS, (B)T:tauS)R:tauR; with A->B introgression at time tauT
 * root age tauR
 * 
 * ## MSCM model ##
 * There are three variants: IM, IIM and SC
 * IIM: (A, (B)T:tauT)R:tauR; with A->T migration during (tauT,tauR)
 * IM:  IIM with tauT = 0
 * SC:  (A, (B)T:tauT)R:tauR; with A->B migration during (0,tauT)
 * 
 */

#include "paml.h"

// CHOOSE direction of inference
// int MtoI = 1;
int MtoI = 0;

// parameters in the migration model (MSCM: IM, IIM, SC)
double w, thetaAm, thetaSm, thetaRm, tauRm, tauTm, thetaTm;

// parameters in the introgression model (MSCI)
double phi, thetaSi, thetaRi, tauSi, tauRi;

// scaling factor for integral tranform: y = exp(-2 / thetaint * t);  t = -thetaint / 2 * y;
double thetaint;

// 0 = IIM or IM, 1 = secondary contact (SC) model
int sc = 0;      

// 0 = IM, 1 = IIM, 2 = SC
int im_model = 0;

// determines which set of parameters to estimate
// MtoI:
//  mode = 0: phi, tauR
//  mode = 1: phi, tauR, tauS
//  mode = 2: phi, thetaS, thetaR
//  mode = 3: phi, thetaS, thetaR, tauR, tauS
//  mode = 4: phi, thetaS = thetaR, tauR, tauS
//
// ItoM:
//  mode = 0: w, tauR
//  mode = 1: w, tauR, tauS (iim and sc only)
//  mode = 2: w, thetaT, thetaR
//  mode = 3: w, thetaT, thetaR, tauR, tauS (iim and sc only)
//  mode = 4: w, thetaT = thetaR, tauR, tauS (iim and sc only)
int mode = 0;

// KL option:
// klt = KL in coalescent time
// klx = KL in number of sites differed between the two sequences
int is_klt = 0;   

// number of sites for klx
int n = 100;

int debug = 0;

int ContourPlot(int mode, int is_klt);

char const* ImModel[] = {"IM", "IIM", "SC"};


// transformation t:(a, infinity) --> y(0, 1)
#define yt(t,a,thetaint)  (1 - exp(-2 / thetaint * (t - a)))
#define ty(y,a,thetaint)  (a - thetaint / 2 * log(1-y))

extern int noisy;


// Transition probabilities of IM model with A->B migration
// Markov chain has states AB, AA, A|B
// Initial state is AB
// Note that if a -> b, then (exp(a) - exp(b))/(a - b) -> exp(a)
void TransitionProbability(double P[], double t, double w, double thetaA) 
{
   double droot = 2 / thetaA - w, ewt = exp(-w * t);
   
   // P00(t)
   P[0] = ewt;

   // P01(t)
   if (fabs(droot) < 1e-20)
   {
      P[1] = w * t * ewt;
   }
   else 
   {
      P[1] = w / droot * (ewt - exp(-2 / thetaA * t));
   }
}


// log density of IIM model in coalescent time (t)
double lnft_m(double t, double w, double thetaT, double thetaR, 
   double tauT, double tauR) 
{
   double lnf=-1e6, Pt[3], Ptau[3];

   TransitionProbability(Pt, t - tauT, w, thetaT);
   TransitionProbability(Ptau, tauR - tauT, w, thetaT);

   if (t > tauT && t < tauR)
   {
      lnf = log(Pt[1] * 2 / thetaT);
   }
   else if (t > tauR && tauR > tauT)
   {
      lnf = log((Ptau[0] + Ptau[1]) * 2 / thetaR) - 2 / thetaR * (t - tauR);
   }

   return(lnf);
}


// log density of SC model in t
double lnft_s(double t, double w, double thetaA, double thetaT, 
   double thetaR, double tauT, double tauR)
{ 
   double lnf=-1e6, Pt[3], Ptau[3], v1, v2, vmax;

   TransitionProbability(Pt, t, w, thetaA);
   TransitionProbability(Ptau, tauT, w, thetaA);

   if (t > 0 && t <= tauT)
   {
      lnf = log(Pt[1] * 2 / thetaA);
   }
   else if (t > tauT && t <= tauR)
   {
      lnf = log(Ptau[1] * 2 / thetaT) - 2 / thetaT * (t - tauT);
   }
   else if (t > tauR && tauR >= tauT)
   {
      v1 = log(Ptau[0]);
      v2 = log(Ptau[1]) - 2 / thetaT * (tauR - tauT);
      vmax = max2(v1, v2);
      lnf = log(exp(v1 - vmax) + exp(v2 - vmax)) + vmax + 
            log(2 / thetaR) - 2 / thetaR * (t - tauR);
   }

   return(lnf);
}


// log density of IIM model in number of sites (x) that differ between two
// sequences assuming a Poisson mutation model
double lnfx_poi_m(int k, int n, double w, double thetaT, double thetaR, 
   double tauT, double tauR)
{
   double lnf, droot = 2 / thetaT - w;
   double nw = 2 * n + w;
   double nT = 2 * n + 2 / thetaT;
   double nR = 2 * n + 2 / thetaR;
   double lnga = LnGamma(k + 1);

   if (fabs(droot) < 1e-20)
   {
      double u1, u2, zz, v1a, v1b, v2a, v2b, v3, vmax;
      u1 = (k+1)*log(tauT) - 2*n*tauT;
      u2 = (k+1)*log(tauR) + w*tauT - nw*tauR;

      v1a = 2*log(w) - log(nw) - lnga + u1;
      v1b = 2*log(w) - log(nw) - lnga + u2;

      zz = log(IncompleteGamma(nw * tauR, k + 1, lnga) - 
               IncompleteGamma(nw * tauT, k + 1, lnga));
      v2a = 2*log(w) - (k+2) * log(nw) + w * tauT + log(k + 1) + zz;
      v2b = 2*log(w) - (k+2) * log(nw) + w * tauT + log(nw * tauT) + zz;

      v3 = -w*(tauR - tauT) + log(1 + w * (tauR - tauT)) - 
           (k+1) * log(nR) + log(2 / thetaR) + 2 / thetaR * tauR +
           lpgamma_upper(nR * tauR, k + 1);

      vmax = max2(max2(max2(v1a, v1b), max2(v2a, v2b)), v3);
      lnf = log(exp(v1a - vmax) - exp(v1b - vmax) +
                exp(v2a - vmax) - exp(v2b - vmax) + 
                exp(v3 - vmax)) + vmax;
 
   }
   else
   {
      double v1, v2, v3, vmax, u1, u2, umax, s;

      v1 = w * tauT - (k + 1) * log(nw) +
           log(IncompleteGamma(nw * tauR, k + 1, lnga) - 
               IncompleteGamma(nw * tauT, k + 1, lnga));
      v2 = 2 / thetaT * tauT - (k + 1) * log(nT) +
           log(IncompleteGamma(nT * tauR, k + 1, lnga) - 
               IncompleteGamma(nT * tauT, k + 1, lnga));
      
      u1 = log(2 / w) - w * (tauR - tauT);
      u2 = log(thetaT) - 2 / thetaT * (tauR - tauT);
      umax = max2(u1, u2);
      s = droot > 0 ? 1 : -1;  // ensure terms are positive before taking log
      
      v3 = log(s * (exp(u1 - umax) - exp(u2 - umax))) + umax +
         2 / thetaR * tauR - (k + 1) * log(nR) - log(thetaR) + 
         lpgamma_upper(nR * tauR, k + 1);
      vmax = max2(max2(v1, v2), v3);
      lnf = log(2*w) - log(s * droot) - log(thetaT) + 
            log(s * exp(v1 - vmax) - s * exp(v2 - vmax) + exp(v3 - vmax)) + vmax;
   }

   lnf += k * log(2 * n);

   return(lnf);
}


// log density of SC model in x assuming a Poisson mutation model
double lnfx_poi_s(int k, int n, double w, double thetaA, double thetaS, 
   double thetaR, double tauT, double tauR)
{
   double lnf, droot = 2 / thetaA - w;
   double nw = 2 * n + w;
   double nA = 2 * n + 2 / thetaA;
   double nR = 2 * n + 2 / thetaR;
   double lnga = LnGamma(k + 1);

   if (fabs(droot) < 1e-20)
   {
      double v1, v2, v3, v4, vmax, u1, u2, umax;
      v1 = -(k+2) * log(nw) + log(k+1) + 
           lpgamma(nw * tauT, k + 1);
      v2 = -nw * tauT + (k+1) * log(tauT) - log(nw) - lnga;
      v3 = log(tauT) - (k+1) * log(nA) + 
           log(IncompleteGamma(nA * tauR, k + 1, lnga) - 
               IncompleteGamma(nA * tauT, k + 1, lnga));

      u1 = -w * tauT;
      u2 = -w * tauR + log(w * tauT);
      umax = max2(u1, u2);
      v4 = log(exp(u1 - umax) + exp(u2 - umax)) + umax - 
        2 * log(w) + log(2/thetaR) +
        2 / thetaR * tauR - (k + 1) * log(nR) +
        lpgamma_upper(nR * tauR, k + 1);
      vmax = max2(max2(v1, v2), max2(v3, v4));
      
      lnf = 2 * log(w) + 
        log(exp(v1 - vmax) - exp(v2 - vmax) + 
            exp(v3 - vmax) + exp(v4 - vmax)) + vmax;
   }
   else 
   {
      double v1, v2, v3, vmax, u1, u2, u3, umax;
      int s = droot > 0 ? 1 : -1;  // sign; ensure terms are positive before taking log

      u1 = lpgamma(nw * tauT, k + 1) - (k+1) * log(nw);
      u2 = lpgamma(nA * tauT, k + 1) - (k+1) * log(nA);
      umax = max2(u1, u2);
      v1 = -log(thetaA) + log(s * exp(u1 - umax) - s * exp(u2 - umax)) + umax;

      u1 = -w * tauT;
      u2 = -2/thetaA * tauT;
      umax = max2(u1, u2);
      v2 = -log(thetaA) + 2/thetaA * tauT - (k+1) * log(nA) +
           log(s * exp(u1 - umax) - s * exp(u2 - umax)) + umax +
           log(IncompleteGamma(nA * tauR, k + 1, lnga) - 
               IncompleteGamma(nA * tauT, k + 1, lnga));

      u1 = -w * tauT - 2 / thetaA * (tauR - tauT);
      u2 = -2 / thetaA * tauT - 2 / thetaA * (tauR - tauT);
      u3 = -w * tauT + log(s * droot) - log(w);
      umax = max2(u1, max2(u2, u3));
      v3 = -log(thetaR) + 2 / thetaR * tauR - (k + 1) * log(nR) +
           log(s * exp(u1 - umax) - s * exp(u2 - umax) + exp(u3 - umax)) + umax +
           lpgamma_upper(nR * tauR, k + 1);
      
      vmax = max2(max2(v1, v2), v3);
      lnf = log(2*w) - log(s * droot) + 
            log(exp(v1 - vmax) + exp(v2 - vmax) + exp(v3 - vmax)) + vmax;
   }

   lnf += k * log(2 * n);

   return(lnf);
}


// log density of MSCI model in t
double lnft_i(double t, double phi, double thetaS, double thetaR, 
   double tauS, double tauR)
{
   double lnf = -1e6;

   if (t > tauS && t < tauR)
   {
      lnf = log(phi * 2 / thetaS) - 2 / thetaS * (t - tauS);
   }
   else if (t > tauR && tauR > tauS)
   {
      lnf = log(1 - phi + phi * exp(-2 / thetaS * (tauR - tauS))) + 
            log(2 / thetaR) - 2 / thetaR * (t - tauR);
   }

   return(lnf);
}


// log density of MSCI model in x assuming a Poisson mutation model
double lnfx_poi_i(int k, int n, double phi, double thetaS, double thetaR, 
   double tauS, double tauR)
{
   double lnf;
   double nS = 2 * n + 2 / thetaS;
   double nR = 2 * n + 2 / thetaR;
   double lnga = LnGamma(k + 1);
   double v1, v2, v3, vmax;

   v1 = log(phi * 2 / thetaS) + 2 / thetaS * tauS - (k+1) * log(nS) +
        log(IncompleteGamma(nS * tauR, k + 1, lnga) - 
            IncompleteGamma(nS * tauS, k + 1, lnga));
   v2 = log(phi * exp(-2 / thetaS * (tauR - tauS)) + 1 - phi) +
        log(2 / thetaR) + 2 / thetaR * tauR - (k+1) * log(nR) +
        lpgamma_upper(nR * tauR, k + 1);
   vmax = max2(v1, v2);
   v3 = log(exp(v1 - vmax) + exp(v2 - vmax)) + vmax;

   lnf = k * log(2 * n) + v3;

   return(lnf);
}


// calculate integrand of KL( p(t|Theta_M) || p(t|Theta_I)) for t in [0,Infinity]
double KLt(double t) 
{
   double m, i, kl;

   m  = sc ? lnft_s(t, w, thetaAm, thetaSm, thetaRm, tauTm, tauRm) : 
             lnft_m(t, w, thetaSm, thetaRm, tauTm, tauRm);
   i  = lnft_i(t, phi, thetaSi, thetaRi, tauSi, tauRi);
   kl = m - i;
   kl = MtoI ? exp(m) * kl : exp(i) * (-kl);

   return(kl);
}


// calculate integrand of KL( p(y|Theta_M) || p(y|Theta_I)) for y in [0, 1]
// y = exp(-2 / thetaint * t)
double KLy(double y)
{
   double t, lnm, lni, kl, out=0;

   if (tauRm > tauTm && tauRi > tauSi)
   {
      t = ty(y, 0, thetaint);
      lnm = sc ? lnft_s(t, w, thetaAm, thetaSm, thetaRm, tauTm, tauRm) : 
                 lnft_m(t, w, thetaSm, thetaRm, tauTm, tauRm);
      lni = lnft_i(t, phi, thetaSi, thetaRi, tauSi, tauRi);
      kl = lnm - lni;
      kl = MtoI ? exp(lnm)*kl : exp(lni)*(-kl);
      if (debug)  printf("  KLy %8.5f %8.5f %8.5f\n", lnm, lni, kl);
      out = kl * thetaint / 2 / (1 - y);
   }

   return(out);
}


// calculate KL using Gauss-Legendre quadrature to integrate over y in [0, 1]
// in five segments (0, t1, t2, t3, t4, infty)
double KLf(double x[], int np)
{
   int npoints = 512, i;
   double y[5 + 1] = { 0, 0, 0, 0, 0, 1 }, kl[5] = {0}, klt, t[5] = {0};
   double tmp;

   if (mode == 0) {
      if (MtoI)
      {
         phi = x[0];
         tauRi = x[1];
      }
      else
      {
         w = x[0];
         tauRm = x[1];
      }
   
   } 
   else if (mode == 1)
   {
      if (MtoI)
      { 
         phi = x[0]; tauRi = x[1]; tauSi = is_klt ? x[2] : x[2] * x[1];
      }
      else
      {   
           w = x[0]; tauRm = x[1]; tauTm = (im_model == 0) ? x[2] : x[2] * x[1];
      }

   } 
   else if (mode == 2)
   {
      if (MtoI)
      { 
         phi = x[0]; thetaSi = x[1]; thetaRi = x[2]; tauRi = x[3];
      }
      else
      { 
           w = x[0]; thetaSm = x[1]; thetaRm = x[2]; tauRm = x[3];
      }
   
   }
   else if (mode == 3)
   {
      if (MtoI)
      { 
         phi = x[0]; thetaSi = x[1]; thetaRi = x[2]; tauRi = x[3]; tauSi = is_klt ? x[4] : x[4] * x[3];
      }
      else
      {   
           w = x[0]; thetaSm = x[1]; thetaRm = x[2]; tauRm = x[3]; tauTm = (im_model == 0) ? x[4] : x[4] * x[3]; 
      }
   
   } 
   else if (mode == 4)
   {
      if (MtoI)
      { 
         phi = x[0]; tauRi = x[1]; tauSi = is_klt ? x[2] : x[2] * x[1]; thetaSi = x[3]; thetaRi = x[3];
      }
      else 
      {   
           w = x[0]; tauRm = x[1]; tauTm = (im_model == 0) ? x[2] : x[2] * x[1]; thetaSm = x[3]; thetaRm = x[3];
      }
   
   }
   else
   {
      error2("invalid mode");
   }
   
   if (debug)
   {
      printf("\n");
      if      (mode == 0)
      { 
         printf("para: %8.5f %8.5f \n", x[0], x[1]); 
      }
      else if (mode == 1)
      { 
         printf("para: %8.5f %8.5f %8.5f \n", x[0], x[1], x[2]);
      }
      else if (mode == 2)
      { 
         printf("para: %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3]);
      }
      else if (mode == 3)
      { 
         printf("para: %8.5f %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3], x[4]); 
      }
      else if (mode == 4)
      { 
         printf("para: %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3]);
      }
   }

   if (sc)
   {
      t[1] = min2(tauTm, tauSi);  
      t[2] = max2(tauTm, tauSi);
      t[3] = min2(tauRm, tauRi);
      t[4] = max2(tauRm, tauRi);

      if (t[3] < t[2])
      {
         tmp = t[2]; t[2] = t[3]; t[3] = tmp;
         tmp = y[2]; y[2] = y[3]; y[3] = tmp;
      }

   }
   else
   {
      t[1] = min2(tauTm, tauSi);  
      t[2] = max2(tauTm, tauSi);
      t[3] = min2(tauRm, tauRi);
      t[4] = max2(tauRm, tauRi);
   }

   y[1] = yt(t[1], 0, thetaint);
   y[2] = yt(t[2], 0, thetaint);
   y[3] = yt(t[3], 0, thetaint);
   y[4] = yt(t[4], 0, thetaint);

   if (debug)  printf("tau: %8.5f %8.5f %8.5f %8.5f\n", tauTm, tauSi, tauRm, tauRi);
   if (debug)  printf("  t: %8.5f %8.5f %8.5f %8.5f\n  y: %8.5f %8.5f %8.5f %8.5f\n", t[1], t[2], t[3], t[4], y[1], y[2], y[3], y[4]);
   
   for (i = 1; i < 5; i++)
   { 
      // can start from i=1 since kl[0] is always 0 for IIM
      if ((i == 0 && fabs(       y[1]) < 1e-20) ||
          (i == 1 && fabs(y[1] - y[2]) < 1e-20) ||
          (i == 2 && fabs(y[2] - y[3]) < 1e-20) ||
          (i == 3 && fabs(y[3] - y[4]) < 1e-20) ||
          (i == 4 && fabs(   1 - y[4]) < 1e-20)) continue;
      kl[i] = NIntegrateGaussLegendre(KLy, y[i], y[i + 1], npoints);
   }
   
   klt = kl[0] + kl[1] + kl[2] + kl[3] + kl[4];

   if (debug)  printf("  kl = %8.8f + %8.8f + %8.8f + %8.8f + %8.8f = %8.8f\n", kl[0], kl[1], kl[2], kl[3], kl[4], klt);
   if (klt<0)  printf("  kl = %8.8f + %8.8f + %8.8f + %8.8f + %8.8f = %8.8f < 0!!\n", kl[0], kl[1], kl[2], kl[3], kl[4], klt);
   
   return klt;
}


// calculate KL( p(x|Theta_M) || p(x|Theta_I)) by summing over x in {0, .., n}
double KLx(double x[], int np) 
{
   double lnm, lni, lnmi, kl=0;

   if (mode == 0) 
   {
      if (MtoI)
      {
         phi = x[0]; tauRi = x[1];
      }
      else
      {
         w = x[0]; tauRm = x[1];
      }
   
   } 
   else if (mode == 1) 
   {
      if (MtoI) 
      { 
         phi = x[0]; tauRi = x[1]; tauSi = x[2] * x[1]; 
      }
      else      
      {   
         w = x[0]; tauRm = x[1]; tauTm = x[2] * x[1]; 
      }

   } 
   else if (mode == 2)
   {
      if (MtoI)
      {
         phi = x[0]; thetaSi = x[1]; thetaRi = x[2]; tauRi = x[3]; 
      }
      else      
      {   
         w = x[0]; thetaSm = x[1]; thetaRm = x[2]; tauRm = x[3]; 
      }
   
   } 
   else if (mode == 3) 
   {
      if (MtoI) 
      { 
         phi = x[0]; thetaSi = x[1]; thetaRi = x[2]; tauRi = x[3]; tauSi = x[4] * x[3]; 
      }
      else
      {
         w = x[0]; thetaSm = x[1]; thetaRm = x[2]; tauRm = x[3]; tauTm = x[4] * x[3]; 
      }
   
   } 
   else if (mode == 4) 
   {
      if (MtoI)
      { 
         phi = x[0]; tauRi = x[1]; tauSi = x[2] * x[1]; thetaSi = x[3]; thetaRi = x[3];
      }
      else 
      {
         w = x[0]; tauRm = x[1]; tauTm = x[2] * x[1]; thetaSm = x[3]; thetaRm = x[3];
      }
   
   } else {
      error2("invalid mode");
   }
   
   if (debug) 
   {
      if      (mode == 0) { printf("\npara: %8.5f %8.5f \n", x[0], x[1]); }
      else if (mode == 1) { printf("\npara: %8.5f %8.5f %8.5f \n", x[0], x[1], x[2]); }
      else if (mode == 2) { printf("\npara: %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3]); }
      else if (mode == 3) { printf("\npara: %8.5f %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3], x[4]); }
      else if (mode == 4) { printf("\npara: %8.5f %8.5f %8.5f %8.5f \n", x[0], x[1], x[2], x[3]); }

      printf(" m: %.3f, %.6f, %.6f, %.6f, %.6f\n",   w, thetaSm, thetaRm, tauTm, tauRm);
      printf(" i: %.6f, %.6f, %.6f, %.6f, %.6f\n", phi, thetaSi, thetaRi, tauSi, tauRi);
      printf(" phi %8.5f, w %8.5f\n", phi, w);
   }

   if (tauRm > tauTm && tauRi > tauSi)
   {
      int dd = ceil((2 * tauRm + thetaRm) * n);  // approx distance between two sequences

      // start from the middle
      // ku counting up to n
      // kd counting down to 0
      int ku = dd;
      int kd = ku - 1;

      while (kd >= 0)
      {
         lnm  = sc ? lnfx_poi_s(kd, n, w, thetaAm, thetaSm, thetaRm, tauTm, tauRm) : 
                     lnfx_poi_m(kd, n, w, thetaSm, thetaRm, tauTm, tauRm);
         lni  = lnfx_poi_i(kd, n, phi, thetaSi, thetaRi, tauSi, tauRi);
         lnmi = MtoI ? exp(lnm) * (lnm - lni) : exp(lni) * (lni - lnm);
         
         if (isnan(lnm) || isnan(lni) || fabs(lnmi) < 1e-20)
         {
            if (debug) printf(" break: %8.6f %8.6f %8.6f\n", lnm, lni, lnmi);
            break;
         }

         kl += lnmi;
         if (debug) printf(" kd: %d, %8.6f %8.6f %8.6f\n", kd, lnm, lni, lnmi);
         kd--;
      }

      while (ku <= n)
      {
         lnm  = sc ? lnfx_poi_s(ku, n, w, thetaAm, thetaSm, thetaRm, tauTm, tauRm) : 
                     lnfx_poi_m(ku, n, w, thetaSm, thetaRm, tauTm, tauRm);
         lni  = lnfx_poi_i(ku, n, phi, thetaSi, thetaRi, tauSi, tauRi);
         lnmi = MtoI ? exp(lnm) * (lnm - lni) : exp(lni) * (lni - lnm);

         if (isnan(lnm) || isnan(lni) || fabs(lnmi) < 1e-20) {
            if (debug) printf(" break: %8.6f %8.6f %8.6f\n", lnm, lni, lnmi);
            break;
         }
         kl += lnmi;
         if (debug) printf(" ku: %d, %8.6f %8.6f %8.6f\n", ku, lnm, lni, lnmi);
         ku++;
      }

      if (debug)  printf("  kl = %8.6f\n", kl);
      
      if (kl < 0)
      {
         printf("  kl = %8.6f < 0!! \n", kl);
      }
   }

   return kl;
}


int ContourPlot(int mode, int is_klt) 
{
   int i, j, np;
   double x[4];
   
   FILE *fout;
   if (is_klt)
   {
      fout = fopen("contour-klt.txt", "w"); 
   }
   else
   { 
      fout = fopen("contour-klx.txt", "w"); 
   }

   printf("true model %s, analysis model %s\n", (MtoI ? "M" : "I"), (MtoI ? "I" : "M"));
   printf("M: w   thetaSm  thetaRm  tauRm\n");
   printf("I: phi thetaSi  thetaRi  tauRi\n");

   thetaint = 0.002;
   w = 0.1 * 4 / 0.01;  thetaSm = 0.002; thetaRm = 0.002;  tauRm = 0.006;
   phi = 0.1;  thetaSi = 0.002; thetaRi = 0.002;  tauRi = 0.006;  tauSi = 0.002;
   n = 1000;
   debug=0;

   if (mode == 0)
   {
      np = 2;
      fprintf(fout, "phi\ttauRi\tKL\n");
      for (phi = 0.005; phi < 0.701; phi+=0.005)
      {
         for (tauRi = 0.001; tauRi < 0.00801; tauRi+=0.0005)
         {
            x[0] = phi;
            x[1] = tauRi;

            if (is_klt)
            { 
               fprintf(fout, "%.6f\t%.6f\t%.6f\n", x[0], x[1], KLf(x, np));
            }
            else
            { 
               fprintf(fout, "%.6f\t%.6f\t%.6f\n", x[0], x[1], KLx(x, np));
            }
         }
      }
   
   }
   else if (mode == 1)
   {
      np = 3;
      fprintf(fout, "phi\ttauRi\ttauSi\tKL\n");
      for (phi = 0.005; phi < 0.701; phi += 0.005)
      {
         for (tauRi = 0.001; tauRi < 0.00801; tauRi += 0.0005)
         {
            for (tauSi = 0.0005; tauSi < tauRi; tauSi += 0.0005)
            {
               x[0] = phi;
               x[1] = tauRi;
               x[2] = tauSi/tauRi;

               if (is_klt)
               { 
                  fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\n", 
                          x[0], x[1], x[1]*x[2], KLf(x, np)); 
               }
               else
               { 
                  fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\n", 
                          x[0], x[1], x[1]*x[2], KLx(x, np)); 
               }
            }
         }
      }

   } 
   else if (mode == 4) 
   {
      np = 4;
      fprintf(fout, "phi\ttauRi\ttauSi\ttheta\tKL\n");
      for (phi = 0.005; phi < 0.901; phi += 0.005)
      {
         for (tauRi = 0.001; tauRi < 0.00801; tauRi += 0.0005)
         {
            for (tauSi = 0.0005; tauSi < tauRi; tauSi += 0.0005)
            {
               for (thetaRi = 0.001; thetaRi < 0.005; thetaRi += 0.0005) 
               {
                  x[0] = phi;
                  x[1] = tauRi;
                  x[2] = tauSi/tauRi;
                  x[3] = thetaRi;
                  thetaSi = thetaRi;

                  if (is_klt)
                  { 
                     fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
                             x[0], x[1], x[1]*x[2], x[3], KLf(x, np)); 
                  }
                  else
                  { 
                     fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", 
                             x[0], x[1], x[1]*x[2], x[3], KLx(x, np)); 
                  }
               }
            }
         }
      }

   } 
   else 
   {
      error2("invalid mode");
   }

   exit(0);
}


int main(int argc, char** argv) 
{

   char* ofname = "out.txt";

   noisy = 3;
   SetSeed(-1, 0);

   printf("true model %s, analysis model %s\n", (MtoI ? "M" : "I"), (MtoI ? "I" : "M"));
   printf("M: w   thetaSm  thetaRm  tauRm\n");
   printf("I: phi thetaSi  thetaRi  tauRi\n");

   // specify mode to run 
   if (argc > 1) 
   {
      char *endp;
      long conv = strtol(argv[1], &endp, 10);

      if (*endp != '\0' || conv > INT_MAX || conv < INT_MIN) 
      {
         error2("invalid mode");
      } 
      else 
      {
         mode = conv;
      }
      if (mode < 0 || mode > 4) error2("invalid mode");
   }
   
   // output filename
   if (argc > 2) 
   {
      ofname = argv[2];
   }

   // klt or klx; if klt, n is not used (Inf)
   // is_klt = 1 means KL in t (infinite sequence length)
   // is_klt = 0 means KL in x (finite sequence length n)
   char* klt_opt = "klx";
   if (argc > 3) 
   {
      klt_opt = argv[3];
      if (strcmp(klt_opt, "klt") == 0) is_klt = 1;
   }

  // choice of IM model: IM, IIM, SC
   if (argc > 4) 
   {
      char *endp;
      long conv = strtol(argv[4], &endp, 10);

      if (*endp != '\0' || conv > INT_MAX || conv < INT_MIN) 
      {
         error2("invalid IM model");
      } 
      else 
      {
         im_model = conv;
      }
      if (im_model < 0 || im_model > 2) error2("invalid IM model");
   }

   // number of sites (n) for klx
   if (argc > 5) 
   {
      char *endp;
      long conv = strtol(argv[5], &endp, 10);

      if (*endp != '\0' || conv > INT_MAX || conv < INT_MIN) 
      {
         error2("invalid n");
      } 
      else 
      {
         n = conv;
      }
   }

   printf("n: %d\n", n);
   printf("output: %s\n", ofname);
   printf("mode: %d\n", mode);
   printf("IM model: %s\n", ImModel[im_model]);
   printf("%s (is_klt: %d)\n", klt_opt, is_klt);

   int i, np, ip;
   double t, t1, t2, x[5] = { 0 }, KL, space[100000];

   // optim bounding box for MtoI
   double xbi[5][2] = { { 0.00001, 0.999}, { 1e-6, 0.5 }, { 1e-6, 0.5 }, { 1e-6, 0.5 }, { 1e-6, 0.5 } };

   // list of M
   double M0[50] = { 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.007, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 1, 2 };
   
   // size of M0
   // int nm = 23;

   // optim bounding box for ItoM
   double xbm[5][2] = { { 1e-6, 5000.}, { 1e-6, 0.5 }, { 1e-6, 0.5 }, { 1e-6, 0.5 }, { 1e-6, 0.5 } };
   
   // list of phi
   double phi0[50] = { 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 };
   
   // size of phi0
   int nm = 11;


   // output file
   FILE *fout = fopen(ofname, "w"), *frub = fopen("rub", "w");

   
   // parameter values for the true model
   thetaTm = 0.01;
   phi = 0.2; thetaSi = 0.002; thetaRi = 0.002; tauRi = 0.004; tauSi = 0.002;
   w = 0.1*4/thetaTm;

   // scaling factor for numerical integration
   thetaint = 0.002;


   if (im_model == 0) 
   {
      // IM model
      if (MtoI)
      { 
         tauTm = 0.0; tauRm = 0.002; tauRi = tauRm; tauSi = tauTm; 
      }
      else
      {
         tauTm = 0.0; tauRm = tauRi;
      }
   
   } 
   else if (im_model == 2) 
   {
      // SC model
      sc = 1;

      if (MtoI) tauSi = 0.0;
   }
   
   // output
   if (mode == 0) 
   {
      // number of parameters to estimate
      np = 2;

      if (MtoI)
      {
         fprintf(fout, "M\tphi\ttauR\tKL\n");
      }
      else
      {
         fprintf(fout, "phi\tw\ttauR\tKL\n");
      }
   
   } 
   else if (mode == 1) 
   {
      np = 3;

      if (MtoI)
      {
         fprintf(fout, "M\tphi\ttauR\ttauS\tKL\n");
      }
      else
      {
         fprintf(fout, "phi\tw\ttauR\ttauT\tKL\n");
      }
   
   } 
   else if (mode == 2) 
   {
      np = 4;

      if (MtoI)
      {
         fprintf(fout, "M\tphi\tthetaS\tthetaR\ttauR\tKL\n");
      }
      else {
         fprintf(fout, "phi\tw\tthetaS\tthetaR\ttauR\tKL\n");
      }

   } 
   else if (mode == 3) 
   {
      np = 5;

      if (MtoI)
      {
         fprintf(fout, "M\tphi\tthetaS\tthetaR\ttauR\ttauS\tKL\n");
      }
      else
      {
         fprintf(fout, "phi\tw\tthetaS\tthetaR\ttauR\ttauT\tKL\n");
      }

   } 
   else if (mode == 4) 
   {
      np = 4;
      
      if (MtoI)
      {
         fprintf(fout, "M\tphi\ttauR\ttauS\ttheta\tKL\n");
      }
      else
      {
         fprintf(fout, "phi\tw\ttauR\ttauT\ttheta\tKL\n");
      }

   }
   else 
   {
      error2("invalid mode");
   }


   // for MtoI, calculate w from M
   double w0[50];
   if (MtoI) for (i = 0; i < nm; i++) w0[i] = M0[i] * 4 / thetaTm;
   
   // initial values: MtoI
   // thetaSi = 0.002; thetaRi = 0.002; tauRi = tauRm; tauSi = sc ? 0 : tauTm;
   // phi = 1 - exp(-4 * M0[0] / thetaTm * (tauRm - tauTm));

   // initial values: ItoM
   thetaTm = 0.01; thetaSm = 0.002; thetaRm = 0.002; tauRm = 0.004; tauTm = 0.0;   // IM
   // thetaTm = 0.01; thetaSm = 0.002; thetaRm = 0.002; tauRm = 0.004; tauTm = 0.002; // IIM
   w = -log(1 - phi) / (tauRi - tauSi);
   thetaAm = thetaSm;  // for SC only

   if (!MtoI && sc) tauTm = 5e-6;

   for (ip = 0; ip < nm; ip++)
   {
      if (MtoI)
      {
         w = w0[ip];
      }
      else
      {
         phi = phi0[ip];
      }

      if (mode == 0) 
      {
         if (MtoI) 
         { 
            x[0] = phi;
            x[1] = tauRi;
            thetaSi = thetaSm;
            thetaRi = thetaRm;
            tauSi = sc ? 0.0 : tauTm; 
            xbi[1][0] = tauSi + 1e-6;
            xbi[1][1] = tauRm * 1.5;
         
         } 
         else 
         { 
            x[0] = w;
            x[1] = tauRm;
            thetaSm = thetaSi;
            thetaRm = thetaRi;
            tauTm = (im_model == 0) ? 0.0 : tauSi;
            xbm[1][0] = tauTm + 1e-6;
            xbm[1][1] = tauRi * 4;
         }
      
      } 
      else if (mode == 1) 
      {
         if (MtoI) 
         { 
            x[0] = phi;
            x[1] = tauRi;
            x[2] = is_klt ? tauSi : tauSi/tauRi;
            thetaSi = thetaSm;
            thetaRi = thetaRm; 
            xbi[1][0] = tauSi + 1e-6;
            xbi[1][1] = tauRm * 1.5;
            xbi[2][0] = is_klt ? (sc ? 1e-7 : (tauTm > 0 ? tauTm : 1e-7)) : 0.001;
            xbi[2][1] = is_klt ? (sc ? 1e-6 : (tauTm > 0 ? tauTm + 1e-6 : 1e-6)) : 0.99;
         
         } 
         else 
         { 
            x[0] = w;

            x[1] = tauRm;
            x[2] = (im_model == 0) ? tauTm : tauTm/tauRm;
            thetaSm = thetaSi;
            thetaRm = thetaRi;
            xbm[1][0] = tauTm + 1e-6;
            xbm[1][1] = tauRi * 4;
            xbm[2][0] = (im_model == 0) ? 1e-7 : 0.001;
            xbm[2][1] = (im_model == 0) ? 1e-6 : 0.99;
         }

      } 
      else if (mode == 2) 
      {
         if (MtoI) 
         { 
            x[0] = phi;
            x[1] = thetaSi;
            x[2] = thetaRi;
            x[3] = tauRi;
            tauSi = sc ? 0 : tauTm;
            xbi[3][0] = tauTm + 1e-6;
            xbi[3][1] = tauRm * 1.5;

         } 
         else 
         { 
            x[0] = w;
            x[1] = thetaSm;
            x[2] = thetaRm;
            x[3] = tauRm;
            tauTm = (im_model == 0) ? 0 : tauSi;
            xbm[3][0] = tauSi + 1e-6;
            xbm[3][1] = tauRi * 2;
         }

      } 
      else if (mode == 3) 
      {
         if (MtoI) 
         { 
            x[0] = phi;
            x[1] = thetaSi;
            x[2] = thetaRi;
            x[3] = tauRi;
            x[4] = is_klt ? tauSi : tauSi/tauRi;
            xbi[3][0] = sc ? 1e-7 : tauTm + 1e-6;
            xbi[3][1] = tauRm * 1.5;
            xbi[4][0] = is_klt ? (sc ? 1e-7 : tauTm + 1e-7) : 0.001;
            xbi[4][1] = is_klt ? (sc ? 1e-6 : tauTm + 1e-6) : 0.99;
         
         } 
         else 
         { 
            x[0] = w;
            x[1] = thetaSm;
            x[2] = thetaRm;
            x[3] = tauRm;
            x[4] = (im_model == 0) ? tauTm : tauTm/tauRm;
            xbm[3][0] = tauSi + 1e-6;
            xbm[3][1] = tauRi * 2;
            xbm[4][0] = (im_model == 0) ? 1e-7 : 0.001;
            xbm[4][1] = (im_model == 0) ? 1e-6 : 0.99;
         }
      
      }
      else if (mode == 4)
      {
         if (MtoI)
         { 
            x[0] = phi;
            x[1] = tauRi;
            x[2] = is_klt ? tauSi : tauSi/tauRi;
            x[3] = thetaRi;
            thetaSi = thetaRi;
            xbi[2][0] = is_klt ? (sc ? 1e-7 : tauTm + 1e-7) : 0.001;
            xbi[2][1] = is_klt ? (sc ? 1e-6 : tauTm + 1e-6) : 0.99;

         }
         else
         { 
            x[0] = w;
            x[1] = tauRm;
            x[2] = (im_model == 0) ? tauTm : tauTm/tauRm;
            x[3] = thetaRm;
            thetaSm = thetaRm;
            xbm[1][0] = tauSi + 1e-6;
            xbm[1][1] = tauRi * 2;
            xbm[2][0] = (im_model == 0) ? 1e-7 : 0.001;
            xbm[2][1] = (im_model == 0) ? 1e-6 : 0.99;
         }

      }
      else
      {
         error2("invalid mode");
      }

      // optional: randomise initial values
      x[0] = 1e-6 + x[0] * (0.5+rndu());
      x[1] = 1e-6 + x[1] * (0.5+rndu());
      x[2] = 1e-6 + x[2] * (0.5+rndu());
      x[3] = 1e-6 + x[3] * (0.5+rndu());
      x[4] = 1e-6 + x[4] * (0.5+rndu());
      
      printf("\ninit & interval:\n%.6f\t%.6f\t%.6f \n%.6f\t%.6f\t%.6f \n%.6f\t%.6f\t%.6f \n%.6f\t%.6f\t%.6f \n%.6f\t%.6f\t%.6f\n", 
         x[0], xbm[0][0], xbm[0][1], 
         x[1], xbm[1][0], xbm[1][1], 
         x[2], xbm[2][0], xbm[2][1], 
         x[3], xbm[3][0], xbm[3][1], 
         x[4], xbm[4][0], xbm[4][1]);

      matout2(stdout, x, 1, np, 9, 6);

      if (is_klt) printf("KL = %9.6f\n", KLf(x, np));
      else        printf("KL = %9.6f\n", KLx(x, np));
      if ((MtoI && w == 0) || !MtoI && phi == 0) break;

      if (MtoI) fprintf(frub, "\n\nM[%d] = %.6f\n", ip, M0[ip]);
      else      fprintf(frub, "\n\nphi[%d] = %.6f\n", ip, phi0[ip]);

      if (is_klt) ming2(frub, &KL, KLf, NULL, x, (MtoI ? xbi : xbm), space, 1e-9, np);
      else        ming2(frub, &KL, KLx, NULL, x, (MtoI ? xbi : xbm), space, 1e-9, np);

      matout2(stdout, x, 1, np, 9, 6);
      printf("min KL = %9.6f\n", KL);

      if (mode == 0)
      {
         fprintf(fout, "%.6f\t%.6f\t%.6f\t%.12f\n", 
            (MtoI ? M0[ip] : phi0[ip]), x[0], x[1], KL);

      }
      else if (mode == 1) 
      {
         fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.12f\n", 
            (MtoI ? M0[ip] : phi0[ip]), x[0], x[1], (MtoI && is_klt) || 
            (!MtoI && im_model == 0) ? x[2] : x[2]*x[1], KL);

      }
      else if (mode == 2) 
      {
         fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.12f\n", 
            (MtoI ? M0[ip] : phi0[ip]), x[0], x[1], x[2], x[3], KL);

      }
      else if (mode == 3)
      {
         fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.12f\n", 
            (MtoI ? M0[ip] : phi0[ip]), x[0], x[1], x[2], x[3], 
            (MtoI && is_klt) || (!MtoI && im_model == 0) ? x[4] : x[4]*x[3], KL);

      }
      else if (mode == 4)
      {
         fprintf(fout, "%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.12f\n", 
            (MtoI ? M0[ip] : phi0[ip]), x[0], x[1], 
            (MtoI && is_klt) || (!MtoI && im_model == 0) ? x[2] : x[2]*x[1], x[3], KL);
      }

      printf("M: w   thetaSm  thetaRm  tauTm  tauRm: %.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",   w, thetaSm, thetaRm, tauTm, tauRm);
      printf("I: phi thetaSi  thetaRi  tauSi  tauRi: %.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", phi, thetaSi, thetaRi, tauSi, tauRi);

      if (MtoI) 
      {
         phi = x[0];
         thetaSi = 0.002;
         thetaRi = 0.002;
         tauRi = tauRm;
         tauSi = sc ? 0 : tauTm;
      }
      else
      {
           w = x[0];
           thetaSm = 0.002;
           thetaRm = 0.002;
           tauRm = tauRi;
           tauTm = im_model == 0 ? 0 : tauSi;
      }

   }
   fclose(fout);
   fclose(frub);
}
