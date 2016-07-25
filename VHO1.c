#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<math.h>
#include<time.h>

#define M 6
#define N 9
#define AP 4

#define Noise_dBm -108 //-80
#define Noise_dBm_AP -86 //-80
#define gammaAP	33	//dBm
#define gammaBS 42	//dBm
#define Power_PIE -114  //dBm
#define Wbs	5000000
#define Wap 22000000
#define MAXLOAD 0.6
#define alpha 0.4

#define size 1000
#define width 100
#define divisions1 2 //x-cord
#define divisions2 1 //y-cord

#define SIR_THS 0.015848 /* -18 dB */
#define K1 4

#define CPICH_RSCP -114

#define PilotPower 0.1
#define Frequency 2000 /* in MHz */ 
#define infinity 1000000.0
#define scale 10

#define HB 10 //hight of base station in meter
#define HM 1 //hight of mobile in meter
#define Gain 15.65
#define gain_AP 20

#define FadeMargin 13	
#define SHOFadeMargin 8
#define SHOTh 3
#define pi 3.141593 

//#define FM 8
//#define onlyAp 13
//#ap antenna gain 6db

FILE *outxx;

int UTx[N]={100,200,300,400,500,600,700,800,900},UTy[N]={100,200,300,400,500,600,700,800,900};
int BTx[M],BTy[M],dg,Deg[M],ATx[AP],ATy[AP];

int Hid[360];

int APassign[AP][N],BSassign[M][N],ABassign[AP][M][N],AAassign[AP][AP][N],BBassign[M][M][N];
int AP_sat[AP],BS_sat[M];

double Horiz[360];
double MaxPowerBS[M],MaxPowerAP[AP];
double W[N][M],S[N][M];
double I_AP[N][AP],Pow_AP[N][AP],I_BS[N][M],Pow_BS[N][M];
double Pow_AA1[N][AP],Pow_AA2[N][AP],Pow_AA[N][AP][AP],I_AA1[N][AP],I_AA2[N][AP];
double Pow_BB1[N][M],Pow_BB2[N][M],Pow_BB[N][M][M],I_BB1[N][M],I_BB2[N][M];
double Pow_AB1[N][AP],Pow_AB2[N][M],Pow_AB[N][AP][M],I_AB1[N][AP],I_AB2[N][M];
double watt,watt1,DB,DB1,PathLoss,Noise,Euclidean,Horizon[N][M];
double SIR_THS1,SIR_THS2,SIR_THS3,SIR_THS4,SIRTHS[N],zeta;

double P_assign_BS[M],R_assign_AP[AP],P_assign_AP[AP];

struct USER
{
	double B[M][5];
	double A[AP][4];//id,distance from user i the power received from that BS or Ap, without and with sho
	double Rreq;
	int degree_A;
	int degree_AA;
	int degree_B;
	int degree_BB;
	int degree_AB;
	int assign;
	int A_sat[AP];
	int B_sat[M];
	int AB_sat[AP][M];
	int AA_sat[AP][AP];
	int BB_sat[M][M];	
}user[N];

void uniform123(int xmin,int xmax,int ymin,int ymax)
{
	int width1, height1, x0, y0, xsep, ysep, i ,j, k, tx;
	
	width1=xmax-xmin;
	height1=ymax-ymin;//?????
	xsep = width1/divisions1;
	ysep=width1/divisions2;
	x0=xsep/2;
	y0=ysep/2;
	
	
	tx=0;
    dg=0;
	for(i=0;i<divisions1;i++){
           for (j=0;j<divisions2;j++){
                            for(k=0;k<3;k++){			
							BTx[tx]=x0+(xsep*i);
                        	BTy[tx]=y0+(ysep*j);
                        	Deg[dg]=120*k;
							tx++;
                        	dg++;
	       	}
	   }
     }
     
    //placement of APs
    x0=xsep/4;
    y0=ysep/4;
	tx=0;
	for(i=0;i<divisions1;i++)
	{
		for(j=0;j<divisions1;j++)
		{
			ATx[tx]=x0+(xsep*j);
			//printf("\n%d",ATx[tx]);
			ATy[tx]=y0+(ysep*i);
			tx++;
		}
	} 
}

void HLoss(void)
{
int i,j;
double epsilon; 


for (i=0;i<N;i++)
{
for (j=0;j<M;j++)
{

if (UTx[i]==BTx[j] && UTy[i]==BTy[j])
{
Horizon[i][j] = Gain;
}
else
{
//First quad
	if((UTx[i]-BTx[j]>=0) && (UTy[i]-BTy[j]>0) && (Deg[j]==0))
	{
		epsilon=90-atan((double)(UTy[i]-BTy[j])/(double)(UTx[i]-BTx[j]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]>=0) && (UTy[i]-BTy[j]>0) && (Deg[j]==120))
	{
		epsilon=330-atan((double)(UTy[i]-BTy[j])/(double)(UTx[i]-BTx[j]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]>=0) && (UTy[i]-BTy[j]>0) && (Deg[j]==240))
	{
		epsilon=210-atan((double)(UTy[i]-BTy[j])/(double)(UTx[i]-BTx[j]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

//Second quad

	if((UTx[i]-BTx[j]<0) && (UTy[i]-BTy[j]>=0) && (Deg[j]==0))
	{
		epsilon=270+atan((double)(UTy[i]-BTy[j])/(double)(BTx[j]-UTx[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]<0) && (UTy[i]-BTy[j]>=0) && (Deg[j]==120))
	{
		epsilon=150+atan((double)(UTy[i]-BTy[j])/(double)(BTx[j]-UTx[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]<0) && (UTy[i]-BTy[j]>=0) && (Deg[j]==240))
	{
		epsilon=30+atan((double)(UTy[i]-BTy[j])/(double)(BTx[j]-UTx[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}


//Third quad

	if((UTx[i]-BTx[j]<=0) && (UTy[i]-BTy[j]<0) && (Deg[j]==0))
	{
		epsilon=180+atan((double)(BTx[j]-UTx[i])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]<=0) && (UTy[i]-BTy[j]<0) && (Deg[j]==120))
	{
		epsilon=60+atan((double)(BTx[j]-UTx[i])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]<=0) && (UTy[i]-BTy[j]<0) && (Deg[j]==240))
	{
		if(atan((double)(BTx[j]-UTx[i])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi) <=60)
		{
			epsilon=300+atan((double)(BTx[j]-UTx[i])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
			//fprintf(outxx,"%d ",(int)epsilon);
			Horizon[i][j]=Horiz[(int)epsilon];
		}
		else
		{
			epsilon=atan((double)(BTx[j]-UTx[i])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi)-60;
			//fprintf(outxx,"%d ",(int)epsilon);
			Horizon[i][j]=Horiz[(int)epsilon];
		}
	}

//Fourth quad

	if((UTx[i]-BTx[j]>0) && (UTy[i]-BTy[j]<=0) && (Deg[j]==0))
	{
		epsilon=180-atan((double)(UTx[i]-BTx[j])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}

	if((UTx[i]-BTx[j]>0) && (UTy[i]-BTy[j]<=0) && (Deg[j]==120))
	{
		if(atan((double)(UTx[i]-BTx[j])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi) <=60)
		{
			epsilon=60-atan((double)(UTx[i]-BTx[j])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
			//fprintf(outxx,"%d ",(int)epsilon);
			Horizon[i][j]=Horiz[(int)epsilon];
		}
		else
		{
			epsilon=420-atan((double)(UTx[i]-BTx[j])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
			//fprintf(outxx,"%d ",(int)epsilon);
			Horizon[i][j]=Horiz[(int)epsilon];
		}
	}

	if((BTx[i]-UTx[j]>0) && (UTy[i]-BTy[j]<=0) && (Deg[j]==240))
	{
		epsilon=300-atan((double)(UTx[i]-BTx[j])/(double)(BTy[j]-UTy[i]))*((double)180/(double)pi);
		//fprintf(outxx,"%d ",(int)epsilon);
		Horizon[i][j]=Horiz[(int)epsilon];
	}
}

}
//fprintf(outxx,"\n\n\n ");
}
/*
i=0;
j=0;
while(i<N)
{
j=0;
while(j<M)
printf("%lf ",Horizon[i][j++]);
printf("\n/>");
i++;

}*/

}

void initiate()
{
	int i,j,k,t;
	double temp;
	//calculate the distances And received powers for each user with and without sho 
	//printf("cordi %d %d %d %d",Tx[1],Ty[1],x[2],y[2]);

	
	for (i=0;i<N;i++)
	{
		user[i].assign=0;
		for(j=0;j<M;j++)
		{
			user[i].B[j][0]=j+1;
			//user[i].B[j][1]= (double)sqrt(pow(x[i]-Tx[j],2)+pow(y[i]-Ty[j],2))/(double)1000; /* Euclidean distance in KM */
			//user[i].B[j][1]=sqrt((Tx[j]-x[i])^2+(Ty[j]-y[i])^2);
			if (UTx[i]==BTx[j] && UTy[i]==BTy[j])
			{
				user[i].B[j][1]=0.02;
				//printf("\n-->%lf",user[i].B[j][1]);
				PathLoss= 69.55+26.16*log10(Frequency)-13.82*log10(HB)-HM*1.1*log10(Frequency)+0.7*HM + 1.56*log10(Frequency)-0.8+(44.9 - 6.55*log10(HB))*log10((double)20/(double)1000);
				//PathLoss=42.6 + 26*log10((double)width/(double)2000) + 20*log10(Frequency);
				user[i].B[j][2]=pow(10,((double)(MaxPowerBS[j] - PathLoss + Gain - Horizon[i][j]- FadeMargin -30)/(double)10));
				user[i].B[j][3]=pow(10,((double)(MaxPowerBS[j] - PathLoss + Gain - Horizon[i][j]- SHOFadeMargin -30)/(double)10));
				user[i].B[j][4]=pow(10,((double)(PathLoss - Gain + Horizon[i][j])/(double)10));		
			}
			else
			{
				//user[i].B[j][1]=0;
				user[i].B[j][1]=(double)sqrt(pow((double)(UTx[i]-BTx[j]),2.0)+pow((double)(UTy[i]-BTy[j]),2.0))/(double)1000; /* Euclidean distance in KM */
				//printf("\n-->%lf",user[i].B[j][1]);
				PathLoss= 69.55+26.16*log10(Frequency)-13.82*log10(HB)-HM*1.1*log10(Frequency)+0.7*HM + 1.56*log10(Frequency)-0.8+(44.9 - 6.55*log10(HB))*log10(user[i].B[j][1]);
				//PathLoss=42.6 + 26*log10(Euclidean) + 20*log10(Frequency);
				user[i].B[j][2]=pow(10,((double)(MaxPowerBS[j] - PathLoss + Gain - Horizon[i][j] - FadeMargin -30)/(double)10));
				user[i].B[j][3]=pow(10,((double)(MaxPowerBS[j] - PathLoss + Gain - Horizon[i][j] - SHOFadeMargin -30)/(double)10));
				user[i].B[j][4]=pow(10,((double)(PathLoss - Gain + Horizon[i][j])/(double)10));	
				//printf("\n---%d",user[i].B[j][1]);
			}
			P_assign_BS[j]=0;
			//printf("\n-->%lf..%d",user[i].B[j][2],Tx[j]);
			//printf("%lf",((double)(MaxPower[j] - PathLoss + Gain - Horizon[i][j] - FadeMargin -30)/(double)10));
			//printf("\nUsers %d(%d,%d) Cij BS%lf(%d,%d):%lf %lf %lf %lf.",i,UTx[i],UTy[i],user[i].B[j][0],BTx[j],BTy[j],user[i].B[j][1],pow(10,scale)*user[i].B[j][2],pow(10,scale)*user[i].B[j][3],user[i].B[j][4]);

		}

		for(j=0;j<AP;j++)
		{
			user[i].A[j][0]=j+1;
			//user[i].A[j][1]=(double)sqrt(pow(x[i]-ATx[j],2)+pow(y[i]-ATy[j],2)); /* Euclidean distance in M */
			//user[i].A[j][1]=sqrt((ATx[j]-x[i])^2+(ATy[j]-y[i])^2);
			if (UTx[i]==ATx[j] && UTy[i]==ATy[j])
			{
				user[i].A[j][1]=20;
				PathLoss= 38+8+28*log10(user[i].A[j][1]);
				//PathLoss=42.6 + 26*log10((double)width/(double)2000) + 20*log10(Frequency);
				user[i].A[j][2]=pow(10,((double)(MaxPowerAP[j] - PathLoss + gain_AP -30)/(double)10));
				//user[i].A[j][3]=pow(10,((double)(MaxPower[j] - PathLoss + Gain - Horizon[i][j] - 30)/(double)10));

			}
			else
			{
				user[i].A[j][1]=(double)sqrt(pow(UTx[i]-ATx[j],2)+pow(UTy[i]-ATy[j],2)); /* Euclidean distance in M */
				PathLoss= 38+8+28*log10(user[i].A[j][1]);
				//PathLoss=42.6 + 26*log10(Euclidean) + 20*log10(Frequency);
				user[i].A[j][2]=pow(10,((double)(MaxPowerAP[j] - PathLoss + gain_AP - FadeMargin -30)/(double)10));
				//user[i].A[j][3]=pow(10,((double)(MaxPower[j] - PathLoss + Gain - Horizon[i][j] - SHOFadeMargin -30)/(double)10));
				//printf("\n---%lf",user[i].A[j][2]);
			}
			P_assign_AP[j]=0;
			//Total_I_AP=+user[i].A[j][2];
			//printf("\nUsers %d(%d,%d) Dij AP%lf(%d,%d):%lf,%lf",i,UTx[i],UTy[i],user[i].A[j][0],ATx[j],ATy[j],user[i].A[j][1],pow(10,scale)*user[i].A[j][2]);
		}

	}
	
	printf("Initiation done!!\n");
}

DegAP(int i)
{
	int j,k;
	user[i].degree_A=0;
	for(j=0;j<AP;j++)
		user[i].A_sat[j]=0;
		
	for(j=0;j<AP;j++)
	{
		//interference cal
		I_AP[i][j]=0;
		for(k=0;k<AP;k++)
		{
			if(k!=j)
			I_AP[i][j]+=P_assign_AP[k]*user[i].A[k][2];
		}	
		Pow_AP[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AP[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/Wap)-1))/user[i].A[j][2]);
		//*pow(10,MaxPowerAP[j]/10)
		if(0.6>=P_assign_AP[j]+Pow_AP[i][j])
		{
		//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			user[i].degree_A++;
			user[i].A_sat[j]=1;
		}
	//	printf("\n user:%d %d",i,user[i].degree_A);
	}
}

DegAP_AP(int i)
{
	int j,k,l;
	double R,R2;
	
	user[i].degree_AA=0;
	for(j=0;j<AP;j++)
	{
		for(k=0;k<AP;k++)
		user[i].AA_sat[j][k]=0;
	}
		
	for(j=0;j<AP;j++)
	{
		//interference cal
		I_AA1[i][j]=0;
		for(k=0;k<AP;k++)
		{
			if(k!=j)
			I_AA1[i][j]+=P_assign_AP[k]*user[i].A[k][2];
		}	
		Pow_AA1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AA1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/Wap)-1))/user[i].A[j][2]);
		//*pow(10,MaxPowerAP[j]/10)
		if(0.6>=P_assign_AP[j]+Pow_AA1[i][j])
		{
		//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			Pow_AA1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AA1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/(Wap*2))-1))/user[i].A[j][2]);
			for(l=j+1;l<AP;l++)
			{
				if(j!=l)
				{
					I_AA2[i][l]=0;
					for(k=0;k<AP;k++)
					{
						if(k!=l)
						I_AA2[i][l]+=P_assign_AP[k]*user[i].A[k][2];
					}	
					Pow_AA2[i][l]=((pow(10,((Noise_dBm_AP-30)/10))+I_AA2[i][l])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/(Wap*2))-1))/user[i].A[l][2]);
					if(0.6>=P_assign_AP[l]+Pow_AA2[i][l])
					{
						user[i].degree_AA++;
						user[i].AA_sat[j][l]=1;
						Pow_AA[i][j][l]=Pow_AA1[i][j]+Pow_AA2[i][l];
					}
				}
			}
		}
		else
		{
			//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			R=Wap*log2(((1*user[i].A[j][2])/(pow(10,(gammaAP-30)/10)*(pow(10,((Noise_dBm_AP-30)/10))+I_AA1[i][j])))+1);//check
			//printf("\n AP: %d Rate:%lf",j,R);
			R2=user[i].Rreq-R;
			Pow_AA1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AA1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,R/Wap)-1))/user[i].A[j][2]);
			for(l=j+1;l<AP;l++)
			{
				if(j!=l)
				{
					I_AA2[i][l]=0;
					for(k=0;k<AP;k++)
					{
						if(k!=l)
						I_AA2[i][l]+=P_assign_AP[k]*user[i].A[k][2];
					}	
					Pow_AA2[i][l]=((pow(10,((Noise_dBm_AP-30)/10))+I_AA2[i][l])*(pow(10,(gammaAP-30)/10)*(pow(2,R2/Wap)-1))/user[i].A[l][2]);
					if(0.6>=P_assign_AP[l]+Pow_AA2[i][l])
					{
						user[i].degree_AA++;
						user[i].AA_sat[j][l]=1;
						Pow_AA[i][j][l]=Pow_AA1[i][j]+Pow_AA2[i][l];
					}
				}
			}
		}
	//	printf("\n AP-AP user:%d %d",i,user[i].degree_AA);
	}
}

DegAP_BS(int i)
{
	int j,k,l;
	double R,R2;
	
	user[i].degree_AB=0;
	for(j=0;j<AP;j++)
	{
		for(k=0;k<M;k++)
		user[i].AB_sat[j][k]=0;
	}
		
	for(j=0;j<AP;j++)
	{
		//interference cal
		I_AB1[i][j]=0;
		for(k=0;k<AP;k++)
		{
			if(k!=j)
			I_AB1[i][j]+=P_assign_AP[k]*user[i].A[k][2];
		}	
		Pow_AB1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AB1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/Wap)-1))/user[i].A[j][2]);
		//*pow(10,MaxPowerAP[j]/10)
		if(0.6>=P_assign_AP[j]+Pow_AB1[i][j])
		{
		//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			Pow_AB1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AB1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,user[i].Rreq/(Wap*2))-1))/user[i].A[j][2]);
			for(l=0;l<M;l++)
			{
					I_AB2[i][l]=0;
					for(k=0;k<M;k++)
					{
						if(k!=l)
						I_AB2[i][l]+=P_assign_BS[k]*user[i].B[k][3];
					}	
					Pow_AB2[i][l]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(pow(10,((Noise_dBm-30)/10))+I_AB2[i][j]))/(double)(user[i].B[j][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(1-alpha))));
					if(0.6>=P_assign_AP[l]+Pow_AB2[i][l])
					{
						user[i].degree_AB++;
						user[i].AB_sat[j][l]=1;
						Pow_AB[i][j][l]=Pow_AB1[i][j]+Pow_AB2[i][l];
					}
			}
		}
		else
		{
			//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			R=Wap*log2(((1*user[i].A[j][2])/(pow(10,(gammaAP-30)/10)*(pow(10,((Noise_dBm_AP-30)/10))+I_AB1[i][j])))+1);//check
			//printf("\n AP: %d Rate:%lf",j,R);
			R2=user[i].Rreq-R;
			Pow_AB1[i][j]=((pow(10,((Noise_dBm_AP-30)/10))+I_AB1[i][j])*(pow(10,(gammaAP-30)/10)*(pow(2,R/Wap)-1))/user[i].A[j][2]);
			for(l=0;l<M;l++)
			{
					I_AB2[i][l]=0;
					for(k=0;k<M;k++)
					{
						if(k!=l)
						I_AB2[i][l]+=P_assign_BS[k]*user[i].B[k][3];
					}	
					Pow_AB2[i][l]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(R2/Wbs))-1)*(pow(10,((Noise_dBm-30)/10))+I_AB2[i][j]))/(double)(user[i].B[j][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(R2/Wbs))-1)*(1-alpha))));
					if(0.6>=P_assign_BS[l]+Pow_AB2[i][l])
					{
						user[i].degree_AB++;
						user[i].AB_sat[j][l]=1;
						Pow_AB[i][j][l]=Pow_AB1[i][j]+Pow_AB2[i][l];
					}
			}
		}
		//printf("\n AP-BS user:%d %d",i,user[i].degree_AB);
	}
}


DegBS(int i)
{
	int j,k;
	user[i].degree_B=0;
	for(j=0;j<M;j++)
		user[i].B_sat[j]=0;
		
	for(j=0;j<M;j++)
	{
		//interference cal
		I_BS[i][j]=0;
		for(k=0;k<M;k++)
		{
			if(k!=j)
			I_BS[i][j]+=P_assign_BS[k]*user[i].B[k][2];
		}	
						
		Pow_BS[i][j]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/Wbs))-1)*(pow(10,((Noise_dBm-30)/10))+I_BS[i][j]))/(double)(user[i].B[j][2]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/Wbs))-1)*(1-alpha))));
		//*pow(10,MaxPowerBS[j]/10)
		if(0.6>=P_assign_BS[j]+Pow_BS[i][j])
		{
		//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerBS[j]/10)*pow(10,10),P_assign_BS[j]*pow(10,10),Pow_BS[i][j]*pow(10,10));
			user[i].degree_B++;
			user[i].B_sat[j]=1;
		}
	}
	//printf("\n user:%d %d",i,user[i].degree_B);
//	user[i].degree_B=0;
}

DegBS_BS(int i)
{
	int j,k,l;
	double R,R2;
	
	user[i].degree_BB=0;
	for(j=0;j<M;j++)
	{
		for(k=0;k<M;k++)
		user[i].BB_sat[j][k]=0;
	}
		
	for(j=0;j<M;j++)
	{
		//interference cal
		I_BB1[i][j]=0;
		for(k=0;k<M;k++)
		{
			if(k!=j)
			I_BB1[i][j]+=P_assign_BS[k]*user[i].B[k][2];
		}	
		Pow_BB1[i][j]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/Wbs))-1)*(pow(10,((Noise_dBm-30)/10))+I_BB1[i][j]))/(double)(user[i].B[j][2]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/Wbs))-1)*(1-alpha))));
		//*pow(10,MaxPowerAP[j]/10)
		if(0.6>=P_assign_BS[j]+Pow_BB1[i][j])
		{
		//printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			I_BB1[i][j]=0;
			for(k=0;k<M;k++)
			{
				if(k!=j)
				I_BB1[i][j]+=P_assign_BS[k]*user[i].B[k][3];
			}	
			Pow_BB1[i][j]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(pow(10,((Noise_dBm-30)/10))+I_BB1[i][j]))/(double)(user[i].B[j][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(1-alpha))));
			for(l=j+1;l<M;l++)
			{
				if(j!=l)
				{
					I_BB2[i][l]=0;
					for(k=0;k<M;k++)
					{
						if(k!=l)
						I_BB2[i][l]+=P_assign_BS[k]*user[i].A[k][3];
					}	
					Pow_BB2[i][l]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(pow(10,((Noise_dBm-30)/10))+I_BB2[i][j]))/(double)(user[i].B[j][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(user[i].Rreq/(Wbs*2)))-1)*(1-alpha))));
					if(0.6>=P_assign_BS[l]+Pow_BB2[i][l] && abs(user[i].B[j][4]-user[i].B[l][4])<=pow(10,SHOTh/10))
					{
						//printf("  -> %d %d %d",i,j,l);
						user[i].degree_BB++;
						user[i].BB_sat[j][l]=1;
						Pow_BB[i][j][l]=Pow_BB1[i][j]+Pow_BB2[i][l];
					}
				//	else
				//	printf("  ->< %d %d %d",i,j,l);
				}
			}
		}
		else
		{
			//	printf("\n%d %d %lf %lf %lf",i,j,0.6*pow(10,MaxPowerAP[j]/10)*pow(10,10),P_assign_AP[j]*pow(10,10),Pow_AP[i][j]*pow(10,10));
			I_BB1[i][j]=0;
			for(k=0;k<M;k++)
			{
				if(k!=j)
				I_BB1[i][j]+=P_assign_BS[k]*user[i].B[k][3];
			}	
			R=Wbs*log2(((1*user[i].B[j][2])/(pow(10,(gammaBS-30)/10)*(pow(10,((Noise_dBm-30)/10))+I_BB1[i][j])))+1);//check
			//printf("\n BS: %d Rate:%lf deg:%d",j,R,user[i].degree_BB);
			R2=user[i].Rreq-R;
			Pow_BB1[i][j]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(R/Wbs))-1)*(pow(10,((Noise_dBm-30)/10))+I_BB1[i][j]))/(double)(user[i].B[j][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(R/Wbs))-1)*(1-alpha))));
			for(l=j+1;l<M;l++)
			{
				if(j!=l)
				{
					I_BB2[i][l]=0;
					for(k=0;k<M;k++)
					{
						if(k!=l)
						I_BB2[i][l]+=P_assign_BS[k]*user[i].B[k][3];
					}	
					Pow_BB2[i][l]=(double)(pow(10,((gammaBS-30)/10))*(pow(2,(R2/Wbs))-1)*(pow(10,((Noise_dBm-30)/10))+I_BB2[i][l]))/(double)(user[i].B[l][3]*(1-(pow(10,((gammaBS-30)/10))*(pow(2,(R2/Wbs))-1)*(1-alpha))));
				//	printf("          %d %d %d %lf",i,j,l,pow(10,10)*abs(user[i].B[j][4]-user[i].B[l][4]));
					
					if(0.6>=P_assign_BS[l]+Pow_BB2[i][l] && abs(user[i].B[j][4]-user[i].B[l][4])<=pow(10,SHOTh/10))
					//if(0.6>=P_assign_BS[l]+Pow_BB2[i][l])
					{
					//	printf("yes");
						user[i].degree_BB++;
						user[i].BB_sat[j][l]=1;
						Pow_BB[i][j][l]=Pow_BB1[i][j]+Pow_BB2[i][l];
					}
				}
			}
		}
	//	printf("\n BS-BS user:%d %d",i,user[i].degree_BB);
	}
}

choice()
{
	int i;
	for(i=0;i<N;i++)
	{
		DegAP(i);
		DegAP_AP(i);
		DegAP_BS(i);
		DegBS(i);
		DegBS_BS(i);
	}
}

int find_min_deg()
{
	int i,min,minval;
	
	minval=user[0].degree_A+user[0].degree_B;
	//first min which is non assigned
	for(i=0;i<N;i++)
	{
		if(user[i].assign==0)
		{
			min=i;
			minval=user[i].degree_A+user[i].degree_B;
			break;
		}
	}
	for(;i<N;i++)
	{
		if(user[i].assign==0 && minval>user[i].degree_A+user[i].degree_B)
		{
			min=i;
			minval=user[i].degree_A+user[i].degree_B;
		}
	} 
	return min;
}

int find_min2_deg()
{
	int i,min,minval;
	
	//minval=user[0].degree_AA+user[0].degree_BB;
	//first min which is non assigned
	for(i=0;i<N;i++)
	{
		if(user[i].assign==0)
		{
			min=i;
			minval=user[i].degree_AA+user[i].degree_BB;
			break;
		}
	}
	for(;i<N;i++)
	{
		if(user[i].assign==0 && minval>user[i].degree_AA+user[i].degree_BB)
		{
			min=i;
			minval=user[i].degree_AA+user[i].degree_BB;
		}
	} 
	return min;
}

int find_min3_deg()
{
	int i,min,minval;
	
	minval=user[0].degree_AB;
	//first min which is non assigned
	for(i=0;i<N;i++)
	{
		if(user[i].assign==0)
		{
			min=i;
			minval=user[i].degree_AB;
			break;
		}
	}
	for(;i<N;i++)
	{
		if(user[i].assign==0 && minval>user[i].degree_AB)
		{
			min=i;
			minval=user[i].degree_AB;
		}
	} 
	return min;
}

int find_min4_deg()
{
	int i,min,minval;
	
	minval=user[0].degree_A+user[0].degree_AA+user[0].degree_B+user[0].degree_BB+user[0].degree_AB;
	//first min which is non assigned
	for(i=0;i<N;i++)
	{
		if(user[i].assign==0)
		{
			min=i;
			minval=user[i].degree_A+user[i].degree_AA+user[i].degree_B+user[i].degree_BB+user[i].degree_AB;
			break;
		}
	}
	for(;i<N;i++)
	{
		if(user[i].assign==0 && minval>user[i].degree_A+user[i].degree_AA+user[i].degree_B+user[i].degree_BB+user[i].degree_AB)
		{
			min=i;
			minval=user[i].degree_A+user[i].degree_AA+user[i].degree_B+user[i].degree_BB+user[i].degree_AB;
		}
	} 
	return min;
}


assign_AP(int min)
{
	int j,minj;
	for(j=0;j<AP;j++)
	{
		if(user[min].A_sat[j]==1)
		{
			minj=j;
			break;
		}
	}
	for(;j<AP;j++)
	{
		if(Pow_AP[min][minj]>=Pow_AP[min][j] && user[min].A_sat[j]==1)
			minj=j;
	}
	P_assign_AP[minj]+=Pow_AP[min][minj];
	AP_sat[minj]=min;
	user[min].assign=1;
	printf(" ->Asat:%d %d AP:%d %d %lf",min,user[min].assign,minj,AP_sat[minj],pow(10,10)*Pow_AP[min][minj]);	
}

assign_BS(int min)
{
	int j,minj;
	for(j=0;j<M;j++)
	{
		if(user[min].B_sat[j]==1)
		{
			minj=j;
			break;
		}
	}
	for(;j<M;j++)
	{
		if(Pow_BS[min][minj]>=Pow_BS[min][j] && user[min].B_sat[j]==1)
			minj=j;
	}
	P_assign_BS[minj]+=Pow_BS[min][minj];
	BS_sat[minj]=min;
	user[min].assign=1;
	printf(" ->BSsat:%d %d BS:%d %d %lf %d",min,user[min].assign,minj,BS_sat[minj],pow(10,10)*Pow_BS[min][minj],user[min].B_sat[minj]);
}

assign_AA(int min)
{
	int j,k,minj,mink,set;
	//printf("%d %d",minj,mink);
	for(j=0,set=0;j<AP && set==0;j++)
	{
		for(k=0;k<AP;k++)
		{
		
			if(user[min].AA_sat[j][k]==1)
			{
				minj=j;
				mink=k;
				set=1;
				break;
			}
		}
	}
	//printf("%d %d",minj,mink);
	for(;j<AP;j++)
	{
		for(;k<AP;k++)
		{
			if(Pow_AA[min][minj][mink]>=Pow_AA[min][j][k] && user[min].AA_sat[j][k]==1)
			{
				minj=j;
				mink=k;
			}
		}
	}
	//printf("%d %d",minj,mink);
	P_assign_AP[minj]+=Pow_AA1[min][minj];
	P_assign_AP[mink]+=Pow_AA2[min][mink];
	user[min].assign=1;
	user[min].AA_sat[minj][mink]=1;
	printf(" ->AAsat:%d %d AP:%d %d %lf %lf",min,user[min].assign,minj,mink,pow(10,10)*Pow_AA1[min][minj],pow(10,10)*Pow_AA2[min][mink]);	
}
	
assign_BB(int min)
{	
	int j,k,minj,mink,set=0;
	
	
	for(j=0;j<M;j++)
	{
		if(set==0)
		{
		
		for(k=0;k<M;k++)
		{
		
			if(user[min].BB_sat[j][k]==1)
			{
				minj=j;
				mink=k;
				set=1;
				break;
			}
		}
		}
		else
			break;
		
	}
	for(;j<M;j++)
	{
		for(;k<M;k++)
		{
			if(Pow_BB[min][minj][mink]>=Pow_BB[min][j][k] && user[min].BB_sat[j][k]==1)
			{
				minj=j;
				mink=k;
			}
		}
	}
	P_assign_BS[minj]+=Pow_BB1[min][minj];
	P_assign_BS[mink]+=Pow_BB2[min][mink];
	user[min].assign=1;
	//user[min].BB_sat[minj][mink]=1;
	printf(" ->BB sat:%d %d BS:%d %d %lf %lf",min,user[min].assign,minj,mink,pow(10,10)*Pow_BB1[min][minj],pow(10,10)*Pow_BB2[min][mink]);
}

assign_AB(int min)
{
	int j,k,minj,mink,set;
	//printf("vvvdfd%d %d",minj,mink);
	
	for(j=0,set=0;j<AP && set==0;j++)
	{
		for(k=0;k<M;k++)
		{
			if(user[min].AB_sat[j][k]==1)
			{
				minj=j;
				mink=k;
				set=1;
				break;
			}
		}
	}
	//printf("%d %d",minj,mink);
	for(;j<AP;j++)
	{
		for(;k<M;k++)
		{
			if(Pow_AB[min][minj][mink]>=Pow_AB[min][j][k] && user[min].AB_sat[j][k]==1)
			{
				minj=j;
				mink=k;
			}
		}
	}
	//printf("%d %d",minj,mink);
	P_assign_AP[minj]+=Pow_AB1[min][minj];
	P_assign_BS[mink]+=Pow_AB2[min][mink];
	user[min].assign=1;
	user[min].AB_sat[minj][mink]=1;
	printf(" ->AB sat:%d %d AP,BS:%d %d %lf %lf",min,user[min].assign,minj,mink,pow(10,10)*Pow_AB1[min][minj],pow(10,10)*Pow_AB2[min][mink]);	
}
assign1()
{
	int i,min,j;
	static int step=0;
	
	//find user with minimum degree
	min=find_min_deg();
	if(user[min].degree_A==0 && user[min].degree_B==0)
		user[min].assign=-1;
	else
	{
		if(user[min].degree_A!=0)
		{
		
		assign_AP(min);
		}
		else
		{
				//printf("yes");
				assign_BS(min);
		}
	}
	
	choice();
	if(step<N-1)
	{
		printf("\nAt STEP %d: Degree A:%d B:%d",step,user[min].degree_A,user[min].degree_B);
		step++;
		assign1();
	}
}

assign2()
{
	int i,min,j;
	static int step=0;
	
	//find user with minimum degree
	min=find_min2_deg();
	//printf("\n%d %d",min,user[min].degree_BB);
	if(user[min].degree_AA==0 && user[min].degree_BB==0)
		user[min].assign=-1;
	else
	{
		if(user[min].degree_AA!=0)
		{
			assign_AA(min);
		}
		else
		{
			printf("***");
			//printf("yes");
			assign_BB(min);
		}
	}
	
	choice();
	if(step<N-1)
	{
		printf("\nAt STEP %d: Degree AA:%d BB:%d",step,user[min].degree_AA,user[min].degree_BB);
		step++;
		assign2();
	}
}

assign3()
{
	int i,min,j;
	static int step=0;
	
	//find user with minimum degree
	min=find_min3_deg();
	printf("\n%d %d",min,user[min].degree_BB);
	if(user[min].degree_AB==0)
		user[min].assign=-1;
	else
	{
		assign_AB(min);
	}
	
	choice();
	if(step<N-1)
	{
		printf("\nAt STEP %d: Degree AB:%d ",step,user[min].degree_AB);
		step++;
		assign3();
	}
}

assign4()
{
	int i,min,j;
	static int step=0;
	
	//find user with minimum degree
	min=find_min_deg();
	printf("\nAt STEP %d: Degree A:%d B:%d AA:%d BB:%d AB:%d",step,user[min].degree_A,user[min].degree_B,user[min].degree_AA,user[min].degree_BB,user[min].degree_AB);
	if(user[min].degree_A==0 && user[min].degree_AA==0 && user[min].degree_B==0 && user[min].degree_BB==0 && user[min].degree_AB==0)
	{
		user[min].assign=-1;
		printf(" %d Not Assigned",min);
	}
	else
	{
		if(user[min].degree_AA!=0 )
		{
			assign_AA(min);
		}
		else
		{
			if(user[min].degree_A!=0 )
			{
				printf("yes");
				assign_AP(min);
			}
			else
			{
				if(user[min].degree_AB!=0)
				{
					printf("yes");
					assign_AB(min);
				}
				else
				{
					if(user[min].degree_BB==0 )
					{
						printf("yes");
						assign_BB(min);
					}
					else
					{
							printf("yes");
							assign_BS(min);
					}
				}
			}
		}
	}
	
	choice();
	if(step<N-1)
	{
		step++;
		assign4();
	}
}

int main(void)
{

FILE *inpx,*inpxx;

int i,j,k,u;
int stp[N];

outxx=fopen("testalpha0.5.lp","w"); 
inpx=fopen("thshold.txt","r");
inpxx=fopen("LD_Horiz.dat","r");
printf("At Main()");

//user data rate assign 
user[0].Rreq=12200000.0;
//user[0].Rreq=12200.0;
user[1].Rreq=164000.0;
//user[1].Rreq=64000.0;
user[2].Rreq=144000.0;
//user[2].Rreq=144000.0;
user[3].Rreq=122000.0;
user[4].Rreq=64000.0;
user[5].Rreq=1440000.0;
user[6].Rreq=12200000000.0;
user[7].Rreq=64000.0;
user[8].Rreq=144000.0; 


for(i=0;i<360;i++)
fscanf(inpxx,"%d %lf\n",&Hid[i],&Horiz[i]);

for(j=0;j<M;j++)
MaxPowerBS[j]=43;

for(j=0;j<AP;j++)
MaxPowerAP[j]=23;

printf("//%lf",0.6*pow(10,MaxPowerAP[0]/10)*pow(10,10));
uniform123(0,size,0,size);
HLoss();
initiate();
choice();
//assign1();	//simple
//assign2(); //with SHO
//assign3(); //with VHO
assign4(); //with simple+VHO+SHO


//location
/*printf("\n");
for(i=0;i<N;i++)
 printf("%d %d:",UTx[i],UTy[i]);

printf("\n");
for(i=0;i<AP;i++)
 printf("%d %d:",ATx[i],ATy[i]); 

printf("\n");
for(i=0;i<M;i++)
 printf("%d %d:",BTx[i],BTy[i]); 
*/

printf("\nUsers alloted:");
for(j=0;j<N;j++)
printf("%d:",user[j].assign);

//check pow assign
for(k=0;k<AP;k++)
{
	printf("\n////AP%lf",pow(10,0)*P_assign_AP[k]);
}
for(k=0;k<M;k++)
{
	printf("\n////BS%lf",pow(10,0)*P_assign_BS[k]);
}	
printf("\nEND");

return(0);
}



