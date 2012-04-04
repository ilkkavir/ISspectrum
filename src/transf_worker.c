#include <stdlib.h>

void Transf(double *pPr,double *nin0Pr,double *tit0Pr,double *mim0Pr,double *psiPr,double *viPr,double *p_m0,int *nion)
{
	register double me,m0,psi0,p;
	register int n;

	/*N/N0*/
	nin0Pr[*nion]=pPr[0];			/* ele */
	for(n=1,p=0.;n<*nion;n++){
		nin0Pr[n]=pPr[0]*pPr[n+4]; p+=nin0Pr[n]; }
	nin0Pr[0]=pPr[0]-p;			/* ion1 */
	/*T/T0*/
	tit0Pr[*nion]=pPr[1]*pPr[2];
	/*m/m0*/
	m0=p_m0[0];
	me=9.1093897e-31/1.6605402e-27;
	mim0Pr[*nion]=me/m0;
	/*v/v0*/
	viPr[*nion]=pPr[4];
	/*psi*/
	psi0=pPr[3]		/*/(2*3.1415926535897932385)*/;
	psiPr[*nion]=psi0*0.35714;
	/*ions*/
	for(n=0;n<*nion;n++){
		tit0Pr[n]=pPr[1];
		mim0Pr[n]=p_m0[n]/m0;
		viPr[n]=pPr[4];
		psiPr[n]=psi0;
	}
}
