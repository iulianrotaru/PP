#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void XORSHIFT32(int numbers,unsigned int *RANDOM,unsigned int R0)
{
	int i;
	unsigned int r = R0;
	RANDOM[0] = r;
    for(i=1;i<numbers;++i)
    {
        r = r ^ (r << 13);
        r = r ^ (r >> 17);
        r = r ^ (r << 5);
        RANDOM[i] = r;
    }
}
unsigned char *liniarizare_poza(char* nume_fisier_sursa,unsigned char *header,unsigned int *dim_img,unsigned int *inaltime_img,unsigned int *latime_img,int *padding)
{
    FILE *fin = fopen(nume_fisier_sursa, "rb");

    if(fin == NULL)
   	{
   		printf("nu am gasit imaginea %s\n",nume_fisier_sursa);
   		return NULL;
   	}

    fseek(fin, 2, SEEK_SET);
    fread(&*dim_img, sizeof(unsigned int), 1, fin);
    ///printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&*latime_img, sizeof(unsigned int), 1, fin);
    fread(&*inaltime_img, sizeof(unsigned int), 1, fin);
    ///printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

    ///calculam padding-ul pentru o linie
	*padding=0;
    if(*latime_img&3)
        *padding=4-(3*(*latime_img))&3;

    fseek(fin,0,SEEK_SET);
    fread(header, 54, 1, fin);
    unsigned char * Q;
    Q = (unsigned char*) calloc(*inaltime_img*(*latime_img)*3,sizeof(unsigned char));
    unsigned int i,j;
    unsigned char pRGB[3];
    for(i = 0; i < *inaltime_img; i++)
	{
		for(j = 0; j < *latime_img; j++)
		{
			fread(pRGB, 3, 1, fin);
			Q[(*inaltime_img-i-1)*(*latime_img)*3+j*3]   = pRGB[0];
            Q[(*inaltime_img-i-1)*(*latime_img)*3+j*3+1] = pRGB[1];
        	Q[(*inaltime_img-i-1)*(*latime_img)*3+j*3+2] = pRGB[2];
		}
		fseek(fin,*padding,SEEK_CUR);
	}
	fclose(fin);
	return Q;
}
void afisare(char *nume_fisier_iesire,unsigned char *Poza,unsigned char *header,unsigned int dim_img,unsigned int inaltime_img,unsigned int latime_img,int padding)
{
    FILE *fout = fopen(nume_fisier_iesire, "wb");

    fwrite(header, 54, 1, fout);
    int i;
    unsigned char RGB[1];
    RGB[0] = 0;
	for(i = 0; i < inaltime_img; i++)
    {
        fwrite(Poza+(inaltime_img-i-1)*latime_img*3,latime_img*3, 1, fout);
        int j;
        for( j = 0 ; j < padding ; ++j)
            fwrite(RGB,1,1,fout);
    }
	fclose(fout);
}
void permut_pixeli_poza(unsigned char *A,unsigned char *B,unsigned int *P,unsigned int inaltime_img,unsigned int latime_img)
{
    int i,j;
    for(i = 0; i < inaltime_img; i++)
		for(j = 0; j < latime_img; j++)
		{
        	/// pozitia in permutare pentru pixelul (i,j)
        	int poz = i*latime_img+j;
        	int ii = P[poz]/latime_img;
        	int jj = P[poz]%latime_img;
        	A[ii*latime_img*3+jj*3]   = B[i*latime_img*3+j*3];
        	A[ii*latime_img*3+jj*3+1] = B[i*latime_img*3+j*3+1];
        	A[ii*latime_img*3+jj*3+2] = B[i*latime_img*3+j*3+2];
		}
}
void read_secret_key(char *nume_cheie,unsigned int *R0,unsigned int *SV)
{
    FILE *fon = fopen(nume_cheie,"r");
	fscanf(fon,"%u %u",&*R0,&*SV);
	fclose(fon);
}
unsigned int *fac_permutarea(char *nume_cheie,unsigned int *SV,unsigned int inaltime_img,unsigned int latime_img,unsigned int Random[])
{
    /// citesc cheia secreta
	unsigned int R0;
	read_secret_key(nume_cheie,&R0,&*SV);

	/// obtin numerele random
    XORSHIFT32(2*inaltime_img*latime_img,Random,R0);

    /// FAC PERMUTAREA
	unsigned int * P;
    P = (unsigned int*) calloc(inaltime_img*latime_img,sizeof(unsigned int));
    int i,j;
	for(i=0;i<inaltime_img*latime_img;++i) P[i]=i;

    int ind = 1;
    for(i=inaltime_img*latime_img-1;i>0;--i)
    {
        int p = Random[ind]%(i+1);
        int aux = P[p];
        P[p] = P[i];
        P[i] = aux;
        ind++;
    }
    return P;
}
void criptare_poza(char *nume_cheie,unsigned char *Poza,unsigned int inaltime_img,unsigned int latime_img)
{
    unsigned int * Random;
    Random = (unsigned int*) calloc(2*inaltime_img*latime_img,sizeof(unsigned int));
	unsigned int SV;
	unsigned int *P = fac_permutarea(nume_cheie,&SV,inaltime_img,latime_img,Random);

	/// permut pixelii
	unsigned char * PozaPermutata;
    PozaPermutata = (unsigned char*) calloc(inaltime_img*latime_img*3,sizeof(unsigned char));
    permut_pixeli_poza(PozaPermutata,Poza,P,inaltime_img,latime_img);

    /// codific poza
    unsigned int act = SV;
    int ind = inaltime_img*latime_img;
    int i,j;
	for(i = 0; i < inaltime_img; i++)
		for(j = 0; j < latime_img; j++)
		{
		    act^=Random[ind];
		    act^=PozaPermutata[i*latime_img*3+j*3];
		    act^=((unsigned int)PozaPermutata[i*latime_img*3+j*3+1]<<8);
		    act^=((unsigned int)PozaPermutata[i*latime_img*3+j*3+2]<<16);
        	Poza[i*latime_img*3+j*3]   =  act     &((1<<8)-1);
        	Poza[i*latime_img*3+j*3+1] = (act>>8) &((1<<8)-1);
        	Poza[i*latime_img*3+j*3+2] = (act>>16)&((1<<8)-1);
        	ind++;
		}

    /// dezaloc memoria
    free(Random);
    free(P);
    free(PozaPermutata);
}
void decriptare_poza(char *nume_cheie,unsigned char *Poza,unsigned int inaltime_img,unsigned int latime_img)
{
    unsigned int * Random;
    Random = (unsigned int*) calloc(2*inaltime_img*latime_img,sizeof(unsigned int));
    unsigned int SV;
	unsigned int *P = fac_permutarea(nume_cheie,&SV,inaltime_img,latime_img,Random);

    /// fac inversa
	unsigned int * invP;
    invP = (unsigned int*) calloc(inaltime_img*latime_img,sizeof(unsigned int));
    int i,j;
    for(i = 0; i < inaltime_img*latime_img ; ++i)
        invP[P[i]]=i;
	/// aplic decodificarea
	unsigned char * PozaSemiDecriptata;
    PozaSemiDecriptata = (unsigned char*) calloc(inaltime_img*latime_img*3,sizeof(unsigned char));

    unsigned int act = SV;
    int ind = inaltime_img*latime_img;
	for(i = 0; i < inaltime_img; i++)
		for(j = 0; j < latime_img; j++)
		{
		    act^=Random[ind];
		    act^=Poza[i*latime_img*3+j*3];
		    act^=((unsigned int)Poza[i*latime_img*3+j*3+1]<<8);
		    act^=((unsigned int)Poza[i*latime_img*3+j*3+2]<<16);
        	PozaSemiDecriptata[i*latime_img*3+j*3]   =  act     &((1<<8)-1);
        	PozaSemiDecriptata[i*latime_img*3+j*3+1] = (act>>8) &((1<<8)-1);
        	PozaSemiDecriptata[i*latime_img*3+j*3+2] = (act>>16)&((1<<8)-1);
        	act = Poza[i*latime_img*3+j*3] + ((unsigned int)Poza[i*latime_img*3+j*3+1]<<8) +((unsigned int)Poza[i*latime_img*3+j*3+2]<<16);
        	ind++;
		}

	/// permut invers pixelii
	permut_pixeli_poza(Poza,PozaSemiDecriptata,invP,inaltime_img,latime_img);

	/// dezaloc memoria
    free(Random);
    free(P);
    free(invP);
    free(PozaSemiDecriptata);
}
void chi_patrat(unsigned char *Poza,unsigned int inaltime_img,unsigned int latime_img)
{
    int ** FR;
    FR = (int**) malloc(3*sizeof(int*));
    int i,j;
    for(i = 0 ; i < 3 ; ++i)
        FR[i] = (int*) calloc(1<<8,sizeof(int));
    for(i = 0; i < inaltime_img; i++)
		for(j = 0; j < latime_img; j++)
        {
            FR[0][Poza[i*latime_img*3+j*3]]++;
            FR[1][Poza[i*latime_img*3+j*3+1]]++;
            FR[2][Poza[i*latime_img*3+j*3+2]]++;
        }
    double valoare_ideala = 1.0 * inaltime_img * latime_img / 256.0;
    double *RGB;
    RGB = (double*) calloc(3,sizeof(double));
    RGB[0]=RGB[1]=RGB[2]=0;
    for( i = 0 ; i < 256 ; ++i)
        for( j = 0 ; j < 3 ; ++j)
            RGB[j]+=1.0*(FR[j][i]-valoare_ideala)*(FR[j][i]-valoare_ideala)/valoare_ideala;
    printf("R : ");
    printf("%.2f\n",RGB[2]);
    printf("G : ");
    printf("%.2f\n",RGB[1]);
    printf("B : ");
    printf("%.2f\n",RGB[0]);

    /// dezaloc memoria
    free(RGB);
    for( i = 0 ; i < 3 ; ++i)
        free(FR[i]);
    free(FR);
}
void criptare_decriptare()
{
    /// citesc nume poza,nume poza criptata de mine, nume imagine criptata, nume imagine decriptata de mine, nume cheie
    char *nume_fisiere = "criptare_decriptare_name_file.txt";
    FILE * fin = fopen(nume_fisiere,"r");

    /// citesc numele pozei
    char * cuvant = NULL;
    char * nume_poza = NULL;
    cuvant = (char *) calloc(50,sizeof(char));
    fscanf(fin,"%s",cuvant);
    int lungime = 0;
    int i;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_poza = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_poza[i] = cuvant[i];
    nume_poza[lungime] = NULL;

    /// citesc numele imaginii pe care o criptez
    char * nume_my_imagine_criptata = NULL;
    fscanf(fin,"%s",cuvant);
    lungime = 0;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_my_imagine_criptata = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_my_imagine_criptata[i] = cuvant[i];
    nume_my_imagine_criptata[lungime] = NULL;

    /// citesc numele imaginii pe care o decriptez
    char * nume_imagine_criptata = NULL;
    fscanf(fin,"%s",cuvant);
    lungime = 0;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_imagine_criptata = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_imagine_criptata[i] = cuvant[i];
    nume_imagine_criptata[lungime] = NULL;

    /// citesc numele imaginii decriptate
    char * nume_imagine_decriptata = NULL;
    fscanf(fin,"%s",cuvant);
    lungime = 0;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_imagine_decriptata = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_imagine_decriptata[i] = cuvant[i];
    nume_imagine_decriptata[lungime] = NULL;

    /// citesc nume cheie
    char * nume_cheie = NULL;
    fscanf(fin,"%s",cuvant);
    lungime = 0;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_cheie = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_cheie[i] = cuvant[i];
    nume_cheie[lungime] = NULL;

    /// inchid fisierul
    fclose(fin);

    /// retin datele pozei
	unsigned int dim_img,latime_img,inaltime_img;
	int padding;
	unsigned char *header;
	header = (unsigned char *) calloc(54,sizeof(unsigned char));
	unsigned char * Poza = liniarizare_poza(nume_poza,header,&dim_img,&inaltime_img,&latime_img,&padding);

	/// afisez pe ecran valorile pentru chi_patrat pentru poza initiala
	printf("Chi-squared test on RGB channels for %s:\n",nume_poza);
	chi_patrat(Poza,inaltime_img,latime_img);

	/// criptez poza
	criptare_poza(nume_cheie,Poza,inaltime_img,latime_img);

	/// afisez poza criptata
	afisare(nume_my_imagine_criptata,Poza,header,dim_img,inaltime_img,latime_img,padding);

	/// afisez pe ecran valorile pentru chi_patrat pentru poza criptata
	printf("Chi-squared test on RGB channels for %s:\n",nume_my_imagine_criptata);
	chi_patrat(Poza,inaltime_img,latime_img);

    unsigned int dim_imgc,latime_imgc,inaltime_imgc;
	int paddingc;
	unsigned char *headerc;
	headerc = (unsigned char *) calloc(54,sizeof(unsigned char));
	unsigned char * Poza_criptata = liniarizare_poza(nume_imagine_criptata,headerc,&dim_imgc,&inaltime_imgc,&latime_imgc,&paddingc);

	/// decriptez poza
	decriptare_poza(nume_cheie,Poza_criptata,inaltime_imgc,latime_imgc);

	/// afisez poza decriptata
	afisare(nume_imagine_decriptata,Poza_criptata,headerc,dim_imgc,inaltime_imgc,latime_imgc,paddingc);

	/// dezaloc memoria
	free(cuvant);
	free(nume_cheie);
	free(nume_fisiere);
	free(nume_imagine_criptata);
	free(nume_imagine_decriptata);
	free(nume_my_imagine_criptata);
	free(nume_poza);
	free(header);
	free(headerc);
	free(Poza);
	free(Poza_criptata);
}
void make_image_grayscale(unsigned char *Poza,unsigned int inaltime_img,unsigned int latime_img)
{
    int i,j;
	for( i = 0 ; i < inaltime_img ; ++i)
        for( j = 0 ; j < latime_img ; ++j)
        {
            unsigned char aux = 0.299*Poza[i*latime_img*3+j*3+2] + 0.587*Poza[i*latime_img*3+j*3+1] + 0.114*Poza[i*latime_img*3+j*3];
            Poza[i*latime_img*3+j*3]   = aux;
            Poza[i*latime_img*3+j*3+1] = aux;
            Poza[i*latime_img*3+j*3+2] = aux;
        }
}
struct fereastra
{
    int x;
    int y;
    int xx;
    int yy;
    int cifra;
    double scor;
};
struct fereastra * template_matching(char * nume_imagine,char * nume_cifra,double prag,int *nr_ferestre,int cf)
{
    unsigned int dim_img,latime_img,inaltime_img;
	int padding_img;
	unsigned char *header_img;
	header_img = (unsigned char *) calloc(54,sizeof(unsigned char));
	unsigned char * Poza = liniarizare_poza(nume_imagine,header_img,&dim_img,&inaltime_img,&latime_img,&padding_img);

	/// transform poza in imagine grayscale
	make_image_grayscale(Poza,inaltime_img,latime_img);

	unsigned int dim_s,latime_s,inaltime_s;
	int padding_s;
	unsigned char *header_s;
	header_s = (unsigned char *) calloc(54,sizeof(unsigned char));
	unsigned char * Sablon = liniarizare_poza(nume_cifra,header_s,&dim_s,&inaltime_s,&latime_s,&padding_s);

	/// transform sablonul in imagine grayscale
	make_image_grayscale(Sablon,inaltime_s,latime_s);

	/// dimensiunea pixelilor sablonului
	int n = inaltime_s * latime_s;
	int i,j;

	/// calculez media frecventelor pixelilor sablonului
	double media_pixeli_sablon = 0;
	long long s = 0;
	for( i = 0 ; i < inaltime_s ; ++i )
        for( j = 0 ; j < latime_s ; ++j)
            media_pixeli_sablon += 1.0 * Sablon[3 * (i*latime_s + j)];
    media_pixeli_sablon /= n;

    /// calculez deviatia pentru sablon
    double deviatie_pixeli_sablon = 0;
    for( i = 0 ; i < inaltime_s ; ++i)
        for( j = 0 ; j < latime_s ; ++j)
            deviatie_pixeli_sablon += 1.0 * (1.0* Sablon[3 * (i*latime_s + j)] - media_pixeli_sablon) * ( 1.0 * Sablon[3 * (i*latime_s + j)] - media_pixeli_sablon);
    deviatie_pixeli_sablon = sqrt( 1.0/(n-1)*deviatie_pixeli_sablon);

    struct fereastra * ans = NULL;
    *nr_ferestre = 0;
    /// trec la gasit ferestre
	for( i = 0 ; i < inaltime_img ; ++i)
        for( j = 0 ; j < latime_img ; ++j)
    {
        /// calculez media frecventelor pixelilor ferestrei
        double media_pixeli_fereastra = 0;
        int k,l;
        for( k = i ; k < i + inaltime_s ; ++k)
            for( l = j ; l < j + latime_s ; ++l)
                if(k < inaltime_img && l < latime_img)
                    media_pixeli_fereastra += 1.0 * Poza[3 * (k*latime_img + l)];
        media_pixeli_fereastra /= n;

        /// calculez deviatia pentru fereastra
        double deviatie_pixeli_fereastra = 0;
        for( k = i ; k < i + inaltime_s ; ++k)
            for( l = j ; l < j + latime_s ; ++l)
                if(k < inaltime_img && l < latime_img)
                    deviatie_pixeli_fereastra += 1.0 * ( 1.0 * Poza[3 * (k*latime_img + l)] - media_pixeli_fereastra) * ( 1.0 * Poza[3 * (k*latime_img + l)] - media_pixeli_fereastra);
                else deviatie_pixeli_fereastra += 1.0 * ( - media_pixeli_fereastra) * (  - media_pixeli_fereastra);
        deviatie_pixeli_fereastra = sqrt(1.0/(n-1)*deviatie_pixeli_fereastra);

        /// calculez corelatia
        double corelatie = 0;
        for( k = i ; k < i + inaltime_s ; ++k)
            for( l = j ; l < j + latime_s ; ++l)
                if(k < inaltime_img && l < latime_img)
                    corelatie += 1.0 * ( 1.0 * Poza[3 * (k*latime_img + l)] - media_pixeli_fereastra) * ( 1.0 * Sablon[3 * ((k-i)*latime_s + (l-j))] - media_pixeli_sablon) / deviatie_pixeli_fereastra / deviatie_pixeli_sablon ;
                else corelatie += 1.0 * ( - media_pixeli_fereastra) * ( 1.0 * Sablon[3 * ((k-i)*latime_s + (l-j))] - media_pixeli_sablon) / deviatie_pixeli_fereastra / deviatie_pixeli_sablon ;
        corelatie /= n;
        if(corelatie >= prag)
        {
            *nr_ferestre += 1 ;
            ans = realloc(ans , *nr_ferestre * sizeof(struct fereastra));
            ans[*nr_ferestre-1].x = i;
            ans[*nr_ferestre-1].y = j;
            ans[*nr_ferestre-1].xx = i + inaltime_s - 1;
            ans[*nr_ferestre-1].yy = j + latime_s - 1;
            ans[*nr_ferestre-1].scor = corelatie;
            ans[*nr_ferestre-1].cifra = cf;
        }
    }

    /// dezaloc memoria
    free(Poza);
    free(Sablon);
    free(header_img);
    free(header_s);

	return ans;
};
int cmp(const void *a,const void *b)
{
    if((((struct fereastra*)a)->scor) < (((struct fereastra*)b)->scor)) return 1;
    if((((struct fereastra*)a)->scor) > (((struct fereastra*)b)->scor)) return -1;
    int A = (((struct fereastra*)a)->scor) * 1e9;
    int B = (((struct fereastra*)b)->scor) * 1e9;
    return B - A;
}
void sortare_ferestre(struct fereastra *Ferestre,int nr_ferestre)
{
    qsort(Ferestre,nr_ferestre,sizeof(struct fereastra),cmp);
}
void desenez_contur_fereastra(unsigned char * Poza,struct fereastra Fereastra,unsigned char *RGB,int inaltime_img,int latime_img)
{
    if(Fereastra.xx > inaltime_img - 1)
        Fereastra.xx = inaltime_img - 1;
    if(Fereastra.yy > latime_img - 1)
        Fereastra.yy = latime_img - 1;
    int j;
    for(j=Fereastra.x;j<=Fereastra.xx;++j)
    {
        int k ;
        for( k = 0 ; k < 3; ++k)
        {
            Poza[3*(j*latime_img+Fereastra.y)+k]  = RGB[k];
            Poza[3*(j*latime_img+Fereastra.yy)+k] = RGB[k];
        }
    }
    for(j=Fereastra.y;j<=Fereastra.yy;++j)
    {
        int k ;
        for( k = 0 ; k < 3; ++k)
        {
            Poza[3*(Fereastra.x*latime_img+j)+k]  = RGB[k];
            Poza[3*(Fereastra.xx*latime_img+j)+k] = RGB[k];
        }
    }
}
void eliminare_non_maxime(struct fereastra *Ferestre, char * ok, int nr_ferestre)
{
    int i,j;
    for( i = 0 ; i < nr_ferestre ; ++i)
        if(ok[i]==0)
    {
        int j;
        int aria_i = (Ferestre[i].xx-Ferestre[i].x+1)*(Ferestre[i].yy-Ferestre[i].y+1);
        for( j = i+1; j < nr_ferestre ; ++j)
            if(ok[j] == 0)
        {
            int aria_intersectiei = 0;
            int x = Ferestre[i].x;
            if(x < Ferestre[j].x)
                x = Ferestre[j].x;
            int y = Ferestre[i].y;
            if(y < Ferestre[j].y)
                y = Ferestre[j].y;
            int xx = Ferestre[i].xx;
            if(xx > Ferestre[j].xx)
                xx = Ferestre[j].xx;
            int yy = Ferestre[i].yy;
            if(yy > Ferestre[j].yy)
                yy = Ferestre[j].yy;
            if(x<=xx&&y<=yy)
                aria_intersectiei = (xx-x+1)*(yy-y+1);
            int aria_j = (Ferestre[j].xx-Ferestre[j].x+1)*(Ferestre[j].yy-Ferestre[j].y+1);
            double suprapunere = (double)aria_intersectiei/(double)(aria_i+aria_j-aria_intersectiei);
            if(suprapunere > 0.2) ok[j] = 1;
        }
    }
}
void pattern_matching()
{
    /// citesc nume poza si sabloane
    char *nume_fisiere = "template_matching_name_file.txt";
    FILE * fin = fopen(nume_fisiere,"r");

    char * cuvant = NULL;
    char * nume_poza = NULL;
    cuvant = (char *) calloc(50,sizeof(char));
    fscanf(fin,"%s",cuvant);
    int lungime = 0;
    int i;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_poza = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_poza[i] = cuvant[i];
    nume_poza[lungime] = NULL;

    char **nume_sablon = NULL;
    nume_sablon = (char**)malloc(10*sizeof(char*));
    for( i = 0 ; i < 10 ; ++i)
    {
        fscanf(fin,"%s",cuvant);
        lungime = 0;
        int j;
        for( j = 0 ; cuvant[j] != NULL ; ++j) lungime = j + 1;
        nume_sablon[i] = (char *)calloc(lungime+1,sizeof(char));
        for( j = 0 ; j < lungime ; ++j)
            nume_sablon[i][j] = cuvant[j];
        nume_sablon[i][lungime] = NULL;
    }
    char *nume_output = NULL;
    fscanf(fin,"%s",cuvant);
    lungime = 0;
    for( i = 0 ; cuvant[i] != NULL ; ++i) lungime = i + 1;
    nume_output = (char *)calloc(lungime+1,sizeof(char));
    for( i = 0 ; i < lungime ; ++i)
        nume_output[i] = cuvant[i];
    nume_output[lungime] = NULL;
    fclose(fin);

    /// adun ferestrele care respecta cerintele
	double prag = 0.5;
	struct fereastra *Ferestre = NULL;
	int nr_ferestre = 0;
	for( i = 0 ; i < 10 ; ++i)
    {
        int add_ferestre;
        struct fereastra *F = template_matching(nume_poza,nume_sablon[i],prag,&add_ferestre,i);
        Ferestre = realloc(Ferestre , (add_ferestre+nr_ferestre)*sizeof(struct fereastra));
        int j;
        for( j = 0 ; j < add_ferestre ; ++j)
            Ferestre[nr_ferestre++] = F[j];
        free(F);
    }

    /// sortez ferestrele descrescator dupa corelatie
    sortare_ferestre(Ferestre,nr_ferestre);

    /// vector pentru eliminat neeliminat
    char *ok;
    ok = (char *) calloc(nr_ferestre,sizeof(char));

    eliminare_non_maxime(Ferestre,ok,nr_ferestre);

    /// culorile pentru chenarul fiecarei cifre
    unsigned char culori[10][3] = {
                                     {  0,  0,255} ,
                                     {  0,255,255} ,
                                     {  0,255,  0} ,
                                     {255,255,  0} ,
                                     {255,  0,255} ,
                                     {255,  0,  0} ,
                                     {192,192,192} ,
                                     {  0,140,255} ,
                                     {128,  0,128} ,
                                     {  0,  0,128}
                                  };

    /// poza pe care desenez detectiile
    unsigned int dim_img,latime_img,inaltime_img;
	int padding;
	unsigned char *header;
	header = (unsigned char *) calloc(54,sizeof(unsigned char));
	unsigned char * Poza = liniarizare_poza(nume_poza,header,&dim_img,&inaltime_img,&latime_img,&padding);

    /// colorarea chenarelor
    int number = 0;
    for( i = 0 ; i < nr_ferestre ; ++i)
        if(ok[i] == 0)
    {
        number ++ ;
        int j;
        unsigned char *RGB = NULL;
        RGB = (unsigned char *) calloc(3,sizeof(unsigned char));
        for( j = 0 ; j < 3 ; ++j)
            RGB[j] = culori[Ferestre[i].cifra][j];
        desenez_contur_fereastra(Poza,Ferestre[i],RGB,inaltime_img,latime_img);
        free(RGB);
    }

    afisare(nume_output,Poza,header,dim_img,inaltime_img,latime_img,padding);

    /// dezaloc memoria
    free(cuvant);
    free(nume_fisiere);
    free(nume_poza);
    free(nume_sablon);
    free(Ferestre);
    free(Poza);
    free(header);
    free(ok);
}
int main()
{
	criptare_decriptare();
    pattern_matching();
	return 0;
}
