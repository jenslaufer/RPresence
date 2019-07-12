#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SingAnal.h"
#ifndef noRPRES
	#include <R_ext/PrtUtil.h>
    #define printf Rprintf
#endif
void abort(FILE *f,char *s) {
    fprintf(f,"\n***************************************************\n");
    fprintf(f,"******************* PRESENCE ERROR ****************\n");
    fprintf(f,"*** %s ***\n",s);
    fprintf(f,"***************************************************\n");
    fprintf(f,"***************************************************\n");
    exit(1);
}

double atof_comma(char *s, int *n) { 
	  char *x;  float v=99.; *n=0;
	  x=strchr(s,','); 
	  if (x!=NULL) x[0]='.'; 
	  *n=sscanf(s,"%f",&v);
	  return((double)v);
}

int chk_real_labels(char *lbl, int *modtype, int i, int dmnum) { char *tmpptr; extern TSingSpec SingSpec;

    tmpptr=strstr(lbl,"_sin"); if (tmpptr != (void*)NULL) { tmpptr[0]='\0';SingSpec.LnkFn[dmnum][i-1]=sinLnk; }
    tmpptr=strstr(lbl,"_lgt"); if (tmpptr != (void*)NULL) { tmpptr[0]='\0';SingSpec.LnkFn[dmnum][i-1]=logitLnk; }
    tmpptr=strstr(lbl,"_mlg"); if (tmpptr != (void*)NULL) { tmpptr[0]='\0';SingSpec.LnkFn[dmnum][i-1]=mlogitLnk; }
    if (strncmp(lbl,"p10(",4)==0 || strncmp(lbl,"p10[",4)==0) SingSpec.FalsePos=1;
    if (strncmp(lbl,"theta'",6)==0 || strncmp(lbl,"th1",3)==0) { *modtype=7;}
    if (strncmp(lbl,"lambda",6)==0 || strncmp(lbl,"lamda",5)==0 || strncmp(lbl,"u",5)==0) { *modtype=4;}
    if ((strncmp(lbl,"delta",5)==0 || strncmp(lbl,"dlta",4)==0) && *modtype<3) { *modtype=5; }
    if (strncmp(lbl,"PsiB",4)==0 || strncmp(lbl,"psiB",4)==0) { *modtype=3; }
    if (*modtype==1 &&  
       (strcmp(lbl,"gam1")==0 || strcmp(lbl,"eps1")==0 || 
        strncmp(lbl,"gam",3)==0 || strncmp(lbl,"eps",3)==0)) { *modtype=2;}
    if (strncmp(lbl,"phi",3)==0 && *modtype==2) { SingSpec.Alt_Param_Checked=1; }
    if (strncmp(lbl,"psiBa",5)==0 && *modtype==3) { SingSpec.Alt_Param_Checked=1; }
    if (strncmp(lbl,"nu",2)==0 && *modtype==3) { SingSpec.Alt_Param_Checked=2; }
    if (strncmp(lbl,"CR0",3)==0) {     *modtype=6; SingSpec.multiseason=true; SingSpec.Alt_Param_Checked=1; }
    if (strncmp(lbl,"psi1-1(",7)==0) { *modtype=6; SingSpec.multiseason=true;}
    if (strncmp(lbl,"psi0(0-1)",9)==0 || strncmp(lbl,"psi01(0)",8)==0) { *modtype=6; SingSpec.Alt_Param_Checked=0; }
    if (strncmp(lbl,"psi0c(0-2)",10)==0) { *modtype=6; SingSpec.Alt_Param_Checked=3; }
    if (strncmp(lbl,"e(",2)==0 || strncmp(lbl,"E(",2)==0 || strncmp(lbl,"Ent(",4)==0) {   *modtype=8; }
    if (strncmp(lbl,"eta",3)==0) {   *modtype=6;  SingSpec.Alt_Param_Checked=2;}
    if (strcmp(lbl,"R")==0 || strncmp(lbl,"th1.01",6)==0) { *modtype=5; }
    if (strncmp(lbl,"eps",3)==0) { SingSpec.multiseason=true;}
	if (strncmp(lbl,"oA",2)==0) { *modtype=11;}
	return(0);		
}

int read_neighbor_file(char *fname) {	extern TSingSpec SingSpec;
	FILE *g=SingSpec.g, *fcov; char s[2048],*t2,*wptr,errstr[1024]; int i,ll,j,jj,sumn;
    fprintf(g,"site neighbor filename:%s\n\nSite Neighbor matrix(1st 30 rows, 60 cols):\n",fname); 
	printf("opening site neighbor file (%s)...\n",fname);
    if ((fcov=fopen(fname,"r"))==NULL) {
        j=sprintf(errstr,"readdesmat: error opening neighbor file(%s)",fname); 
        abort(g,errstr);
    }		
    for (i=0; (t2=fgets(s,sizeof(s),fcov)); i++) {  
		if (SingSpec.Verbose>0 || i<30) fprintf(g,"%3d:",i); 
        ll=SingSpec.N; if (strlen(s)<(unsigned)SingSpec.N) ll=strlen(s);
        if (i<SingSpec.N) {
            wptr=strchr(s,' '); SingSpec.neighbor_wgt[i]=1;
            if (wptr != (void*)NULL) { ll=wptr-s; SingSpec.neighbor_wgt[i]=atof_comma(wptr,&j); }
            for (j=jj=sumn=0; j<ll; j++) {               //  must be csv or tab-delim file... 
                if (s[j]=='0' || s[j]=='1') {
                    if (SingSpec.Verbose>0 || (i<30 && j<60)) fprintf(g,"%c",s[j]);
                    if (jj<SingSpec.N) {
						SingSpec.neighbor[i][jj]=s[j]-'0';  //  skip every other character
						sumn+=SingSpec.neighbor[i][jj++];
					}
                }
            }
            if (sumn>0) SingSpec.neighbor_wgt[i]/=(double)sumn;
            if (jj<SingSpec.N) SingSpec.neighbor[i][jj]=2; 
            if (SingSpec.Verbose>0 || i<30) fprintf(g," %f\n",SingSpec.neighbor_wgt[i]);
        }
    }  //  end of site i loop
    fclose(fcov); printf("done reading neighbor file\n");	
	return(0);
}

int check_desmat_rowlabels(int *modtype) { extern TSingSpec SingSpec;
    FILE *g=SingSpec.g; int i,ii,j,k,rowflg,dm_has_eps=0,dm_has_gam=0,dm_num_psi=0; char tstr[1024],*tmpptr;
    for (i=ii=0; i<6; i++)                           //  check each row of design matrix
        for (j=0; j<SingSpec.NrowsDM[i]; j++,ii++) {     
            strcpy(tstr,SingSpec.Realname[i][j]); 	
			if (strstr(tstr,"eps")!=NULL) dm_has_eps=1;
			if (strstr(tstr,"gam")!=NULL) dm_has_gam=1;
			if (strstr(tstr,"psi")!=NULL) dm_num_psi++;
            if (strstr(tstr,"link")!=NULL) {
                if (strstr(tstr,"sin")!=NULL) SingSpec.LnkFn[i][j]=sinLnk;
                if (strstr(tstr,"loglog")!=NULL) SingSpec.LnkFn[i][j]=loglogLnk;
                if (strstr(tstr,"exp")!=NULL) SingSpec.LnkFn[i][j]=expLnk;
                if (strstr(tstr,"id")!=NULL) SingSpec.LnkFn[i][j]=IDLnk;
                if (strstr(tstr,"sin/2")!=NULL) SingSpec.LnkFn[i][j]=sinLnk2;
            }
            if ((tmpptr=strchr(SingSpec.Realname[i][j],'='))!=NULL)  {
                SingSpec.fixed[ii]=atof(tmpptr+1);
                fprintf(g,"*%s parm(%d) fixed to %f.\n",SingSpec.Realname[i][j],ii,SingSpec.fixed[ii]);
            } 
            for (k=rowflg=0; k<SingSpec.NParKK[i]; k++)     
                if (SingSpec.DMat[i][j][k]!=0. || SingSpec.DMat_ptr[i][j][k]!=0) rowflg=1;   //  if all cols of row are zero
            if (rowflg==0 && SingSpec.fixed[ii] < -998) {   //  and parm is not already fixed
				SingSpec.fixed[ii]=0.;                      //  then fix it to zero.
                fprintf(g,"*%s parm(%d) fixed to 0.\n",SingSpec.Realname[i][j],ii+1);
            }
        }
    if (*modtype==2 && SingSpec.NrowsDM[0]>1 && SingSpec.Model==1) {
        if (dm_has_eps==0) SingSpec.Model=2;  //  no eps in design matrix -> psi and gamma
        if (dm_has_gam==0) SingSpec.Model=3;  //  no gam in design matrix -> psi and eps
        if (SingSpec.Model==2 && dm_num_psi==1) SingSpec.Model=4;  //  no eps in DM and only 1 psi -> eps=1-gam
    }
	return(0);
}

void readdesmat(int *NPar, double **Params, FILE *f, int *modtype) {
    int i,j,k=0,l,ii=0,jj=0,irows=2,icols=2,dmnum=0,nest=0,nbeta=0,*desmat,*desmatrow,j2,lastmodtype=0;
    char tstr[2560],errstr[256],*tstr2,lbl[256],**bname,*t2,sdm[1024];
	const char *autocov[]={"psi1","psiA","upsi","psiB"};
    double xdm,xval=0; FILE *g;  int get_modl_num(char *s);
    extern TSingSpec SingSpec; g=SingSpec.g; 
    
    *NPar=0; t2=fgets(tstr,sizeof(tstr),f);
    
    while (!feof(f)) {
        if (strstr(tstr,"model=") != (void*)NULL) { 
            *modtype=get_modl_num(tstr); 
			fprintf(g,"\n==>%s\n",tstr); 
			t2=fgets(tstr,sizeof(tstr),f);
        }
        if (tstr[0]<'0' || tstr[0]>'6') break;
        for (i=0; i<((int)strlen(tstr)); i++) if (tstr[i]==',') tstr[i]=' ';  // convert "," to " " in string 
        j=sscanf(tstr,"%d %d %d",&i,&irows,&icols);
        fprintf(g,"\nMatrix %d: rows=%d, cols=%d\n",dmnum+1,irows,icols);
        if (SingSpec.Verbose>1) printf("Matrix %d: rows=%d, cols=%d\n",dmnum+1,irows,icols);
        if (dmnum>5) {
            l=sprintf(errstr,"readdesmat: dmnum=%d tstr=%s i=%d irows=%d icols=%d",dmnum,tstr,i,irows,icols); 
            abort(g,errstr);
        }
        SingSpec.NParKK[dmnum]=0;
        SingSpec.DMat[dmnum] = new double*[irows+1];
        SingSpec.DMat_ptr[dmnum] = new int*[irows+1];
        SingSpec.Realname[dmnum]=new char*[irows+1]; 
        SingSpec.LnkFn[dmnum]=new int[irows+1];
        for (i=0; i<=irows; i++) {
            SingSpec.DMat[dmnum][i]=new double[icols+1]; SingSpec.DMat_ptr[dmnum][i]=new int[icols+1]; 
            for (j=0; j<=icols; j++) SingSpec.DMat[dmnum][i][j]=SingSpec.DMat_ptr[dmnum][i][j]=0;
            SingSpec.Realname[dmnum][i]=new char[64]; strcpy(SingSpec.Realname[dmnum][i],""); 
            SingSpec.LnkFn[dmnum][i]=0;
        }
        SingSpec.Betaname[dmnum]=new char*[icols+1]; 
        SingSpec.BetaFixed[dmnum]=new double[icols+1];
        bname = new char*[icols+1]; 
        for (i=0; i<=icols; i++) {
            SingSpec.Betaname[dmnum][i]=new char[1024]; 
            SingSpec.Betaname[dmnum][i][0]='\0';
            SingSpec.BetaFixed[dmnum][i]=1.1e44; 
            bname[i] = new char[512]; strcpy(bname[i],"");
        }
        if(irows>0 && !feof(f)) {                //   get betanames (if start w/ 'z', then fix them
            t2=fgets(tstr,sizeof(tstr),f); fprintf(g,"            %s",tstr);
            tstr2=strtok(tstr,","); 
            if (SingSpec.Verbose>2) printf("tstr=%s tstr2=%s\n",tstr,tstr2);
            for (j=0; tstr2!=NULL; j++) {
                if (SingSpec.Verbose>2) printf("j=%d tstr2=%s\n",j,tstr2);
                if (j>0 && j<=icols && (strlen(tstr2)>1 || tstr2[0]!='-')) {
                    if (tstr2[0]>='-' && tstr2[0]<='9') {
                        if (j>=icols) {
                            l=sprintf(errstr,"readdesmat: j=%d icols=%d tstr=%s tstr2=%s",j,icols,tstr,tstr2);
                            abort(g,errstr);
                        }
                        if (SingSpec.Verbose>0) printf("dmnum=%d j=%d tstr2=%s\n",dmnum,j,tstr2);
                        SingSpec.BetaFixed[dmnum][j-1]=atof_comma(tstr2,&k);
                    }
                }
                tstr2=strtok(NULL,",");
            }
        }
        for (i=1; i<irows; i++) { 
            if (!feof(f)) {
                t2=fgets(tstr,sizeof(tstr),f); lastmodtype=*modtype;
                if (SingSpec.Verbose>2) printf("tstr:%s\n",tstr);
                j=strchr(tstr,'\n')-tstr; tstr[j]='\0'; tstr2=strtok(tstr,","); fprintf(g,"%-12s ",tstr2);
                if (strlen(tstr2)>=sizeof(lbl)) tstr2[sizeof(lbl)-1]='\0';
                strcpy(lbl,tstr2); 
				j=chk_real_labels(lbl, modtype, i, dmnum);
                strcpy(SingSpec.Realname[dmnum][i-1],lbl);  // save label
				j=chk_real_labels(lbl,modtype, i, dmnum); 
                if (SingSpec.Verbose>1) printf("realname(%d,%d)=%s  icols=%d modtype=%d\n",dmnum,i-1,lbl,icols,*modtype);
				if (*modtype != lastmodtype) printf("realname(%d,%d)=%s  icols=%d modtype=%d\n",dmnum,i-1,lbl,icols,*modtype);
            }
            for (j=1; j<icols; j++) { 
                if (!feof(f)) tstr2=strtok(NULL,","); 
                else {  printf("eof in readdesmat\n"); strcpy(tstr,"1"); tstr2=&tstr[0]; }
                if (tstr2==NULL) { strcpy(tstr,"0"); tstr2=strtok(tstr,",");  }
       
                if (i>=irows || j>=icols) {
                    l=sprintf(errstr,"readdesmat: Dmat(%d %d %d)=%s",dmnum,i-1,j-1,tstr2); 
                    abort(g,errstr);
                }
                xdm=atof_comma(tstr2,&k);   // k=number of valid numbers read from string (=1 if number, =0 if chars)
			    if (k==1) { 
				    SingSpec.DMat[dmnum][i-1][j-1]=xdm; SingSpec.DMat_ptr[dmnum][i-1][j-1]=0;
					if (xdm==0 || xdm==1) sprintf(sdm,"%1.0f",xdm); else sprintf(sdm,"%f",xdm);
				}
				else {
                  //   check if site/srvy covariate
                  for (k=0; k<SingSpec.NSiteCov; k++)  if (strcmp(tstr2,SingSpec.CovNames[k])==0) break;//  check list of site covars
                  if (k<SingSpec.NSiteCov) { SingSpec.DMat_ptr[dmnum][i-1][j-1]=k+1; strcpy(sdm,SingSpec.CovNames[k]);}
				  else {
                    for (k=0; k<SingSpec.NSampCov; k++) if (strcmp(tstr2,SingSpec.CovNames2[k])==0) break;//  check list of sample covars
                    if (k<SingSpec.NSampCov) {	SingSpec.DMat_ptr[dmnum][i-1][j-1]=-(k+1); strcpy(sdm,SingSpec.CovNames2[k]);}
					else {
                        for (k=0; k<4; k++)  if (strcmp(tstr2,autocov[k])==0) break;
						if (k<4) {    // *autocovs[4]={"psi1","psiA","upsi","psiB"};
				    	    SingSpec.Model=5; strcpy(sdm,autocov[k]); 
							SingSpec.DMat_ptr[dmnum][i-1][j-1]=3001+(k>2);
							SingSpec.UseNeighbor=2;
							if (k==2) SingSpec.uncondpsi=1;
							//if (k<2) SingSpec.UseNeighbor=2;
						}
					    else { 
						    printf("ERROR: covar name not found (%s)\n",tstr2); fprintf(g,"ERROR: covar name not found (%s)\n",tstr2);
                            exit(1);
						}
                    }
                  }
				}
                if (SingSpec.DMat_ptr[dmnum][i-1][j-1]==0) {//  if design matrix entry (row i, col j) is not a covariate...
					//if (xdm==0 || xdm==1) fprintf(g," %1.0f",xdm); else fprintf(g," %f",xdm);
				    if (xdm!=0) sprintf(bname[j],"%s.%c%d",SingSpec.Realname[dmnum][i-1],'a'+dmnum,j);
				}
				else  {    
                    if(strlen(bname[j])<1) { //  if betaname is null... assign 1st real parm name to betaname
					   if (strncmp(sdm,SingSpec.Realname[dmnum][i-1],strlen(SingSpec.Realname[dmnum][i-1]))!=0) {
						   strcpy(bname[j],SingSpec.Realname[dmnum][i-1]); strcat(bname[j],"."); 
					   }
					   strcat(bname[j],sdm);
					}
				}
				fprintf(g," %s",sdm);
            }
            fprintf(g,"\n");
        }
        SingSpec.NrowsDM[dmnum]=irows-1; if (irows<1) SingSpec.NrowsDM[dmnum]=0;
        fprintf(g,"========================\n");
        for (i=0; i<(icols-1); i++) {
            if (SingSpec.BetaFixed[dmnum][i]>1.0e44) (*NPar)++;
            else fprintf(g,"Beta(%d,%d) fixed to %f\n",dmnum,i,SingSpec.BetaFixed[dmnum][i]);
			if (strlen(bname[i+1])<1) sprintf(bname[i+1],"%c%d",'a'+dmnum,i+1);
            strncpy(SingSpec.Betaname[dmnum][i],bname[i+1],32); 
        }
        for (i=0; i<=icols; i++) delete[] bname[i]; 
        delete [] bname;
        nest+=SingSpec.NrowsDM[dmnum]; if(icols>0) SingSpec.NParKK[dmnum]=icols-1; 
        t2=fgets(tstr,sizeof(tstr),f);
		if (SingSpec.Verbose>1) {
			for (i=0; i<SingSpec.NrowsDM[dmnum]; i++) {
				fprintf(g,"%-13s",SingSpec.Realname[dmnum][i]);
				for (j=0; j<SingSpec.NParKK[dmnum]; j++) fprintf(g," %1.1f(%d)",SingSpec.DMat[dmnum][i][j],SingSpec.DMat_ptr[dmnum][i][j]);
				fprintf(g,"\n");
			}
		}
        dmnum++;
    }  //  end of loop for design matrix, dmnum
	if (*modtype==7 && SingSpec.NrowsDM[4]>0) *modtype=10;
    nbeta=*(NPar); desmat=new int[nest]; desmatrow=new int[nest]; // desmat,desmatrow used to get design-matrix number and 
    for (i=jj=j2=0; i<6; i++) {                            // row within DM from absolute real parameter number
        for (j=0; j<SingSpec.NrowsDM[i]; j++) { desmat[jj]=i; desmatrow[jj++]=j;} 
    }
    
    SingSpec.fixed = new double[nest]; for (i=0; i<nest; i++) SingSpec.fixed[i]=-999;
    double *Par; Par = new double[nbeta+1]; for (i=0; i<=nbeta; i++) Par[i]= 0.0;
    *Params=Par;

    //   read fixed and init values
    l=strlen(tstr);
    while (l>0) { 
        tstr[l-1]='\0';
        if (l>1) {
            if (strchr(tstr,'(') == NULL) { //    *** get site neighbors from separate file if specified ****
                if (strncasecmp(tstr,"neighbor",8) == 0) j=read_neighbor_file(strchr(tstr,'-')+1);
            }
            else {
                i=atoi(strstr(tstr,"(")+1); ii=desmat[i-1]; jj=desmatrow[i-1];
                if (strchr(tstr,'=')!=NULL) {
                    if (strstr(tstr,")=eq")!=NULL) xval=10.0; //  flag for theta to be computed as equalib
                    else xval=atof_comma(strstr(tstr,"=")+1,&j);
                }
                if (strstr(tstr,"fix")!=NULL) {
                    SingSpec.fixed[i-1]=xval; 
                    fprintf(g,"fix(%d)=%f /* %s */\n",i,xval,SingSpec.Realname[ii][jj]);
                    printf("fix(%d)=%f /* %s */\n",i,xval,SingSpec.Realname[ii][jj]);
                }
                if (strstr(tstr,"init")!=NULL) {
                    if (i<1 || i>=nest)
                        fprintf(g,"\nreaddesmat: init par(1>%d or %d>%d) tstr=%s",i-1,i-1,nest,tstr);
					else {
                        if (strstr(tstr,"real")!=NULL) {
                            if (xval<.0000000001) xval=.0000000001; 
                            if (xval>.9999999999) xval=.9999999999; 
                            Par[i-1]=log(xval/(1.-xval)); fprintf(g,"init(%d)=%f\n",i,Par[i-1]); 
						}
                        else Par[i-1]=xval;
                    }
                }
            }  //  end else
        }   //  end if l>1
        if ((t2=fgets(tstr,sizeof(tstr),f))==NULL) break; l=strlen(tstr);
    }  // end while
    j=check_desmat_rowlabels(modtype);
    if (SingSpec.Verbose>0) printf("done with readdesmat, model=%d modtype=%d\n",SingSpec.Model,*modtype); 
	delete [] desmat; delete [] desmatrow;
}
