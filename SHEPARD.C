/* Gerador de tons de Shepard */
/* By Antonio Carlos M. de Queiroz - acmq@coe.ufrj.br */
/* Version 1.0 - 6/2/97 */

#include <string.h>
#include <math.h>
#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include "xview.h"

#define s_int short int
#define l_int long int

#define maxf 100

/* ".wav" files have a header of 44 bytes: */

typedef struct header {
  char riff[4];             	/* "RIFF" */
  l_int filesize;		/* size of the file */
  char wave[4];              	/* "WAVE" */
  char fmt[4];               	/* "fmt " */
  l_int fmtsize;                /* 16, in general */
  s_int wFormatTag;		/* 1 for PCM */
  s_int nChannels;		/* 1 for mono, 2 for stereo */
  l_int nSamplesPerSec;	        /* sampling frequency (Hz) */
  l_int nAvgBytesPerSec;	/* = nBlockAlign*nSamplesPerSec */
  s_int nBlockAlign;		/* = wBitsPerSample/8 * nChannels */
  s_int wBitsPerSample;		/* 16 or 8 */
  char data[4];		        /* "data" */
  l_int datasize;		/* size of the data block (bytes) */
} header;

  FILE *saida;
  l_int i,j,amostras;
  double dpi,ln10s20,amplitude,gmax,gmin,v,w,t,wmin,wmax,wreal,lnr,nreal;
  s_int vi,n[maxf],sobe;
  header h;
  char buf[100];

/* Declaration of the interface objects */
Xv_opaque frame1,ttaxa,tfmin,tn,tperiodo,
  ttempo,trazao,tarquivo,bexec,message1;

#define set_sel(xx,yy) xx->v.ssetting.sel_setting=yy
#define set_int(xx,yy) xx->v.stextfield.panel_int=yy
#define set_real(xx,yy) xx->v.stextfield.panel_real=(double)yy
#define set_value(xx,yy) strcpy(xx->v.stextfield.panel_value,yy)
#define set_label(xx,yy) xv_set(xx,yy)
#define set_max(xx,yy) xx->v.stextfield.max_value=yy
#define get_dx(xx) xx->dx
#define get_dy(xx) xx->dy
#define get_real(xx) (double)xx->v.stextfield.panel_real
#define get_int(xx) xx->v.stextfield.panel_int
#define get_value(xx) xx->v.stextfield.panel_value
#define get_sel(xx) xx->v.ssetting.sel_setting
#define update(xx) xv_set(xx,xx->xv_label)

/* Callback procedures */

void Acionar(Xv_opaque obj)
{
  amplitude=30000/get_int(tn);
  gmax=0;
  gmin=-60;
  xv_set(message1,"Escrevendo...");
  dpi=2.0*M_PI;
  ln10s20=log(10)/20;
  amostras=(l_int)(get_real(ttempo)*get_real(ttaxa));
  strcpy(h.riff,"RIFF");
  h.filesize=sizeof(header)-8+amostras*2;
  strcpy(h.wave,"WAVE");
  strcpy(h.fmt,"fmt ");
  h.fmtsize=16;
  h.wFormatTag=1;
  h.nChannels=1;
  h.nSamplesPerSec=(int)get_real(ttaxa);
  h.nAvgBytesPerSec=(int)get_real(ttaxa)*16;
  h.nBlockAlign=2;
  h.wBitsPerSample=16;
  strcpy(h.data,"data");
  h.datasize=amostras*2;
  saida=fopen(get_value(tarquivo),"wb");
  fwrite(&h,sizeof(header),1,saida);
  lnr=log(get_real(trazao));
  wmin=dpi*get_real(tfmin);
  wmax=dpi*get_real(tfmin)*exp(get_int(tn)*lnr);
  sobe=wmax>wmin;
  for (i=1; i<=get_int(tn); i++)
    n[i]=i-1;
  for (i=1; i<=amostras; i++) {
    v=0;
    t=i/get_real(ttaxa);
    for (j=1; j<=get_int(tn); j++) {
      nreal=t/get_real(tperiodo)+n[j];
      wreal=wmin*exp(lnr*nreal);
      if ((sobe && wreal>wmax) || (!sobe && wreal<wmax)) n[j]-=get_int(tn);
      w=wmin/t/lnr*get_real(tperiodo)*exp(lnr*(t/get_real(tperiodo)+n[j]));
      /* v+=amplitude*exp(ln10s20*(gmin+(gmax-gmin)*0.5*(1-cos(nreal/get_int(tn)*dpi))))*cos(w*t);*/
      /* v+=amplitude*0.5*(1-cos(nreal/get_int(tn)*dpi))*cos(w*t); */
      v+=amplitude*cos(w*t);
    }
    vi=(s_int)v;
    fwrite(&vi,2,1,saida);
  }
  fclose(saida);
  sprintf(buf,"Arquivo %s escrito.",get_value(tarquivo));
  xv_set(message1,buf);
}

void lista_fmax(Xv_opaque obj)
{
  double fmax;
  fmax=get_real(tfmin)*exp(get_int(tn)*log(get_real(trazao)));
  sprintf(buf,"fmax=%g",fmax);
  if (fmax>get_real(ttaxa)/2) message1->fore_color=RED;
  xv_set(message1,buf);
  message1->fore_color=BLACK;
}

void main()
{
  /* Inicialization */
  xv_init(0,0);
  normal_length=10;
  frame1=xv_create(frame);
    strcpy(frame1->xv_label,"Gerador de tons de Shepard");
    frame1->x=159;
    frame1->y=158;
    frame1->dx=319;
    frame1->dy=159;
  ttaxa=xv_create(textfield);
    strcpy(ttaxa->xv_label,"Taxa de amostragem (Hz)");
    ttaxa->v.stextfield.field_type=real_field;
    ttaxa->v.stextfield.panel_real=11025;
  tfmin=xv_create(textfield);
    strcpy(tfmin->xv_label,"Frequencia minima (Hz)");
    tfmin->y=15;
    tfmin->v.stextfield.field_type=real_field;
    tfmin->v.stextfield.panel_real=100;
    tfmin->notify_handler=lista_fmax;
  tn=xv_create(textfield);
    strcpy(tn->xv_label,"Numero de harmonicos");
    tn->y=30;
    tn->v.stextfield.field_type=int_field;
    tn->v.stextfield.panel_int=4;
    tn->notify_handler=lista_fmax;
  tperiodo=xv_create(textfield);
    strcpy(tperiodo->xv_label,"Periodo da variacao (s)");
    tperiodo->y=45;
    tperiodo->v.stextfield.field_type=real_field;
    tperiodo->v.stextfield.panel_real=1;
  ttempo=xv_create(textfield);
    strcpy(ttempo->xv_label,"Tempo total (s)");
    ttempo->y=60;
    ttempo->v.stextfield.field_type=real_field;
    ttempo->v.stextfield.panel_real=10;
  trazao=xv_create(textfield);
    strcpy(trazao->xv_label,"Razao entre frequencias");
    trazao->v.stextfield.field_type=real_field;
    trazao->v.stextfield.panel_real=2;
    trazao->y=75;
    trazao->notify_handler=lista_fmax;
  tarquivo=xv_create(textfield);
    strcpy(tarquivo->xv_label,"Arquivo de saida");
    tarquivo->y=90;
    strcpy(tarquivo->v.stextfield.panel_value,"tones.wav");
  bexec=xv_create(button);
    strcpy(bexec->xv_label,"Acionar");
    bexec->y=105;
    bexec->notify_handler=Acionar;
  message1=xv_create(message);
    message1->y=120;
  /* Opens all the windows */
  open_window(frame1);
  xv_main_loop(frame1);
  /* Exit */
  restorecrtmode();
}
