/*************************************
 An ANSI-C implementation of the digamma-function for real arguments based
 on the Chebyshev expansion proposed in appendix E of 
 http://arXiv.org/abs/math.CA/0403344 . This is identical to the implementation
 by Jet Wimp, Math. Comp. vol 15 no 74 (1961) pp 174 (see Table 1).
 For other implementations see
 the GSL implementation for Psi(Digamma) in
 http://www.gnu.org/software/gsl/manual/html_node/Psi-_0028Digamma_0029-Function.html
 
 Richard J. Mathar, 2005-11-24
 **************************************/

//#include "BLCA_required_libs.h"
#include "BLCA_digamma.h"

// edited version of original file by Richard J. Mathar for BayesLCA package: edited by J. Wyse 

/** The digamma function in long double precision.
 * @param x the real value of the argument
 * @return the value of the digamma (psi) function at that point
 * @author Richard J. Mathar
 * @since 2005-11-24
 */
long double BLCA_digammal(long double x)
{
  /* force into the interval 1..3 */
  if( x < 0.0L )
    return BLCA_digammal(1.0L-x)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	/* reflection formula */
  else if( x < 1.0L )
    return BLCA_digammal(1.0L+x)-1.0L/x ;
  else if ( x == 1.0L)
    return -M_GAMMAl ;
  else if ( x == 2.0L)
    return 1.0L-M_GAMMAl ;
  else if ( x == 3.0L)
    return 1.5L-M_GAMMAl ;
  else if ( x > 3.0L)
    /* duplication formula */
    return 0.5L*(BLCA_digammal(x/2.0L)+BLCA_digammal((x+1.0L)/2.0L))+M_LN2l ;
  else
  {
    /* Just for your information, the following lines contain
     * the Maple source code to re-generate the table that is
     * eventually becoming the Kncoe[] array below
     * interface(prettyprint=0) :
     * Digits := 63 :
     * r := 0 : 
     * 
     * for l from 1 to 60 do
     * 	d := binomial(-1/2,l) :
     * 	r := r+d*(-1)^l*(Zeta(2*l+1) -1) ;
     * 	evalf(r) ;
     * 	print(%,evalf(1+Psi(1)-r)) ;
     *o d :
     * 
     * for N from 1 to 28 do
     * 	r := 0 :
     * 	n := N-1 :
     *
     *	for l from iquo(n+3,2) to 70 do
     *		d := 0 :
     *		for s from 0 to n+1 do
     *		 d := d+(-1)^s*binomial(n+1,s)*binomial((s-1)/2,l) :
     *		od :
     *		if 2*l-n > 1 then
     *		r := r+d*(-1)^l*(Zeta(2*l-n) -1) :
     *		fi :
     *	od :
     *	print(evalf((-1)^n*2*r)) ;
     *od :
     *quit :
     */
    /*  static long double Kncoe[] = { .30459198558715155634315638246624251L,
                                       .72037977439182833573548891941219706L, -.12454959243861367729528855995001087L,
                                       .27769457331927827002810119567456810e-1L, -.67762371439822456447373550186163070e-2L,
                                       .17238755142247705209823876688592170e-2L, -.44817699064252933515310345718960928e-3L,
                                       .11793660000155572716272710617753373e-3L, -.31253894280980134452125172274246963e-4L,
                                       .83173997012173283398932708991137488e-5L, -.22191427643780045431149221890172210e-5L,
                                       .59302266729329346291029599913617915e-6L, -.15863051191470655433559920279603632e-6L,
                                       .42459203983193603241777510648681429e-7L, -.11369129616951114238848106591780146e-7L,
                                       .304502217295931698401459168423403510e-8L, -.81568455080753152802915013641723686e-9L,
                                       .21852324749975455125936715817306383e-9L, -.58546491441689515680751900276454407e-10L,
                                       .15686348450871204869813586459513648e-10L, -.42029496273143231373796179302482033e-11L,
                                       .11261435719264907097227520956710754e-11L, -.30174353636860279765375177200637590e-12L,
                                       .80850955256389526647406571868193768e-13L, -.21663779809421233144009565199997351e-13L,
                                       .58047634271339391495076374966835526e-14L, -.15553767189204733561108869588173845e-14L,
                                       .41676108598040807753707828039353330e-15L, -.11167065064221317094734023242188463e-15L } ;*/
    
    long double Tn_1 = 1.0L ;	/* T_{n-1}(x), started at n=1 */
    long double Tn = x-2.0L ;	/* T_{n}(x) , started at n=1 */
    long double resul = Kncoe[0] + Kncoe[1]*Tn ;
    
    x -= 2.0L ;
    
    int n;
    for(n = 2 ; n < sizeof(Kncoe)/sizeof(long double) ;n++)
    {
      const long double Tn1 = 2.0L * x * Tn - Tn_1 ;	/* Chebyshev recursion, Eq. 22.7.4 Abramowitz-Stegun */
    resul += Kncoe[n]*Tn1 ;
    Tn_1 = Tn ;
    Tn = Tn1 ;
    }
    return resul ;
  }
}

/** The optional interface to CREASO's IDL is added if someone has defined
 * the cpp macro export_IDL_REF, which typically has been done by including the
 * files stdio.h and idl_export.h before this one here.
 */
//#ifdef export_IDL_REF
/** CALL_EXTERNAL interface.
 * A template of calling this C function from IDL  is
 * @verbatim
 * dg = CALL_EXTERNAL('digamma.so',X)
 * @endverbatim
 * @param argc the number of arguments. This is supposed to be 1 and not 
 *    checked further because that might have negative impact on performance.
 * @param argv the parameter list. The first element is the parameter x
 *    supposed to be of type DOUBLE in IDL 
 * @return the return value, again of IDL-type DOUBLE
 * @since 2007-01-16
 * @author Richard J. Mathar
 */
/*double digamma_idl(int argc, void *argv[])
 {
 long double x = *(double*)argv[0] ;
 return (double)digammal(x) ;
 }
#endif*/ /* export_IDL_REF */
 
 //#ifdef TEST
 
 /* an alternate implementation for test purposes, using formula 6.3.16 of Abramowitz/Stegun with the 
  first n terms */
 //#include <stdio.h>
 
 /*long double digammalAlt(long double x, int n)
  {
  // force into the interval 1..3 
  if( x < 0.0L )
  return digammalAlt(1.0L-x,n)+M_PIl/tanl(M_PIl*(1.0L-x)) ;	// reflection formula 
  else if( x < 1.0L )
  return digammalAlt(1.0L+x,n)-1.0L/x ;
  else if ( x == 1.0L)
  return -M_GAMMAl ;
  else if ( x == 2.0L)
  return 1.0L-M_GAMMAl ;
  else if ( x == 3.0L)
  return 1.5L-M_GAMMAl ;
  else if ( x > 3.0L)
  return digammalAlt(x-1.0L,n)+1.0L/(x-1.0L) ;
  else
  {
  x -= 1.0L ;
  register long double resul = -M_GAMMAl ;
  
  for( ; n >= 1 ;n--)
  resul += x/(n*(n+x)) ;
  return resul ;
  }
  }*/
 
 /*int main(int argc, char *argv[])
  {
  for( long double x=0.01 ; x < 5. ; x += 0.02)
  printf("%.2Lf %.30Lf %.30Lf %.30Lf\n",x, digammal(x), digammalAlt(x,100), digammalAlt(x,200) ) ;
  }*/
 
 //#endif /* TEST */
 