#include "SingAnal.h"
#include <math.h>
#include <stdlib.h>

double MS1Like(double *Params, int NPar)
{
    double getParam(int param_num, int param_grp, int row, int site, int srvy, double *Params, int transf, int LinkFn);
    int ii, jj, seen;
    double psi, p, psi2, p2, dlta, sum, like=0,like1=0,like2=0;
    extern TSingSpec SingSpec;
/*      * computes estimates under multi-strata PRESENCE model 
     *               with 2 strata (nonbreeders=1; breeders=2)
     *      psi1(i) =Pr(patch i is occupied)
     *      psi2(i) =Pr(patch i is occupied by breeders|occupied)
     *      p1(ij)  =Pr(detection in period j | patch i is occupied by nonbreeders)
     *      p2(ij)  =Pr(detection in period j | patch i is occupied by breeders)
     *      dlta(ij)  =Pr(identified as breeders in period j | patch i is occupied
     *                                 by breeders and animals detected)
     * 
     *                                 N
     *                        (1-psi1)/ \  (psi1)
     *                   #not occupied   #occupied
     *                  /                /    \
     *                 /       (1-psi2) /      \  (psi2)
     *                /      #non-breeders    #breeders         
     *               /        /\                / \
     *              /  (1-p1)/  \(p1)    (1-p2)/   \  (p2) 
     *             /     #not #detected     #not   #detected
     *            /    detected    |      detected        /\
     *           /        |        |         |           /  \
     *          /         |        |         |  (1-dlta)no   \ (dlta)  
     *         /          |        |         |      yng-det  yng-det
     *
     *        0           0        1         0          1      2
     *
     */
    // go through data and evaluate likelihood
    for (ii=0, sum=0.; ii<SingSpec.N; ii++) {
        psi=getParam(0,0,0,ii,0,Params,0,SingSpec.LinkFn);
        psi2=getParam(1,0,1,ii,0,Params,0,SingSpec.LinkFn);
        like1=1; like2=1; seen=0;
        for (jj=0; jj<SingSpec.T; jj++) {
            p=getParam(jj+2,1,jj,ii,jj,Params,0,SingSpec.LinkFn);
            p2=getParam(jj+2+SingSpec.T,1,jj+SingSpec.T,ii,jj,Params,0,SingSpec.LinkFn);
            dlta=getParam(jj+2+2*SingSpec.T,2,jj,ii,jj,Params,0,SingSpec.LinkFn);
            if (SingSpec.Data[ii][jj]!=-1) {
                if (SingSpec.Data[ii][jj]==1) { 
                    like1 *= p; 
                    like2 *= p2*(1-dlta); 
                    if (seen==0) seen = 1;
                } 
                if (SingSpec.Data[ii][jj]==2) { 
                    like2 *= p2*dlta; 
                    if (seen<2) seen = 2; 
                } 
                if (SingSpec.Data[ii][jj]==0) {
                    like1 *= (1.0-p);        
                    like2 *= (1.0-p2);        
                }
            }
        }
        like=psi*psi2*like2;
        if (seen<2) like += psi*(1.-psi2)*like1;
        if (seen==0) like += (1.0-psi);
        if (like<SingSpec.nearzero) like=SingSpec.nearzero;
        sum -= SingSpec.det_hist_frq[ii]*log(like); 
		SingSpec.expval[ii]=like;
        //printf("%d %d%d%d%d%d%d %f %f %f %f %f %f\n",ii+1,SingSpec.Data[ii][0],SingSpec.Data[ii][1],SingSpec.Data[ii][2],
        //       SingSpec.Data[ii][3],SingSpec.Data[ii][4],SingSpec.Data[ii][5],psi,psi2,p,p2,like,sum);
    }
    // delete dynamic variables
    //printf("sum=%19.14f",sum); for (ii=0;ii<9;ii++) printf(" %12.7g",Params[ii]); printf("\n");
    return(sum);
}
