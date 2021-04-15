/* Waveform Maker - Generates several signals as wave files */
/* By Antonio Carlos M. de Queiroz - acmq@coe.ufrj.br */
/* Version 1.0 - 1/12/96 */
/* Version 1.1 - 25/1/97 - Works with BC */

#include <string.h>
#include <math.h>
#include <values.h>
#include "xview.h"
#include "xv_var.h"

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

#define s_int short int
#define l_int long int

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
  l_int i,amostras,bloco;
  double dpi,a1,b1,a2,b2,a3,b3,wa1,wa2,wa3,wt1,wt2,ph1,ph2,t,amp1,amp2,amp3,amp41,amp42;
  s_int v1,v2,v3,v4,vn,v,an;
  header h;

/* Declaration of the interface objects */
Xv_opaque fmain,tstart1,tstart2,tstart3,tend1,
  tend2,tend3,tampl1,tampl2,tampl3,
  tampln,tname,tsampling,tsamples,bwrite,
  ttrw1,ttrw2,ttramp1,ttramp2,ttrfase1,ttrfase2,tbloco,
  message1,bquit;

/* Callback procedures */

void NotifyWrite(Xv_opaque obj)
{
  /* Notify handler for bwrite */
  char buf[100];

  dpi=2.0*M_PI;
  amostras=(l_int)get_real(tsamples);
  strcpy(h.riff,"RIFF");
  h.filesize=sizeof(header)-8+(l_int)get_real(tsamples)*2;
  strcpy(h.wave,"WAVE");
  strcpy(h.fmt,"fmt ");
  h.fmtsize=16;
  h.wFormatTag=1;
  h.nChannels=1;
  h.nSamplesPerSec=(l_int)get_real(tsampling);
  h.nAvgBytesPerSec=(l_int)get_real(tsampling)*16;
  h.nBlockAlign=2;
  h.wBitsPerSample=16;
  strcpy(h.data,"data");
  h.datasize=(l_int)get_real(tsamples)*2;
  saida=fopen(get_value(tname),"wb");
  fwrite(&h,sizeof(header),1,saida);
  a1=dpi*(get_real(tend1)-get_real(tstart1))/get_real(tsamples);
  b1=dpi*get_real(tstart1)-a1;
  a2=dpi*(get_real(tend2)-get_real(tstart2))/get_real(tsamples);
  b2=dpi*get_real(tstart2)-a2;
  a3=dpi*(get_real(tend3)-get_real(tstart3))/get_real(tsamples);
  b3=dpi*get_real(tstart3)-a3;
  amp1=get_int(tampl1);
  amp2=get_int(tampl2);
  amp3=get_int(tampl3);
  an=get_int(tampln);
  wt1=dpi*get_real(ttrw1);
  wt2=dpi*get_real(ttrw2);
  ph1=M_PI/180*get_real(ttrfase1);
  ph2=M_PI/180*get_real(ttrfase2);
  bloco=(l_int)get_real(tbloco);
  amp41=get_int(ttramp1);
  amp42=get_int(ttramp2);
  for (i=1; i<=amostras; i++) {
    wa1=a1/2*i+b1; /* /2! */
    wa2=a2/2*i+b2;
    wa3=a3/2*i+b3;
    t=i/get_real(tsampling);
    v1=(int)(amp1*cos(wa1*t));
    v2=(int)(amp2*cos(wa2*t));
    v3=(int)(amp3*cos(wa3*t));
#ifdef __GNUC__
    vn=2*(int)((double)random()/MAXINT*an)-an;
#else
    vn=2*random(an)-an;
#endif
    if ((i/bloco)&1) v4=(int)(amp42*cos(wt2*t+ph2));
    else v4=(int)(amp41*cos(wt1*t+ph1));
    v=v1+v2+v3+vn+v4;
    fwrite(&v,2,1,saida);
  }
  fclose(saida);
  sprintf(buf,"File %s written.",get_value(tname));
  xv_set(message1,buf);
}

void NotifyQuit(Xv_opaque obj)
{
  /* Notify handler for bquit */
  xv_end=1;
}

void main()
{
  /* Inicialization */
  xv_init(0,0);
  /* Menu creation */
  /* Interface objects creation */
  fmain=xv_create(NULL,FRAME,
    XV_LABEL,"Wave Maker",
    X,19,
    Y,10,
    DX,599,
    DY,189,
    NULL);
  tstart1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Sweep 1: Start",
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    NULL);
  tstart2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Sweep 2: Start",
    Y,15,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"5000",
    NULL);
  tstart3=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Sweep 3: Start",
    Y,30,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"4000",
    NULL);
  tend1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"End",
    X,205,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"5000",
    NULL);
  tend2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"End",
    X,205,
    Y,15,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"1000",
    NULL);
  tend3=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"End",
    X,205,
    Y,30,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    NULL);
  tampl1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Amplitude",
    X,322,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,5000,
    MIN_VALUE,0,
    NULL);
  tampl2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Amplitude",
    X,322,
    Y,15,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,5000,
    MIN_VALUE,0,
    NULL);
  tampl3=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Amplitude",
    X,322,
    Y,30,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,5000,
    MIN_VALUE,0,
    NULL);
  tampln=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Noise: Amplitude",
    Y,45,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,1000,
    MIN_VALUE,0,
    NULL);
  ttrw1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Transition: Freq. 1",
    Y,60,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"500",
    NULL);
  ttrw2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Freq. 2",
    Y,60,
    X,250,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"1000",
    NULL);
  ttramp1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Ampl. 1",
    Y,75,
    X,96,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,5000,
    MIN_VALUE,0,
    NULL);
  ttramp2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Ampl. 2",
    Y,75,
    X,250,
    VALUE_LENGTH,10,
    FIELD_TYPE,int_field,
    PANEL_INT,5000,
    MIN_VALUE,0,
    NULL);
  ttrfase1=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Phase 1",
    Y,90,
    X,96,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NULL);
  ttrfase2=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Phase 2",
    Y,90,
    X,250,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NULL);
  tbloco=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Sa./block",
    Y,90,
    X,400,
    VALUE_LENGTH,10,
    FIELD_TYPE,real_field,
    PANEL_REAL,"10000",
    NULL);
  tname=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"File name (.wav)",
    Y,120,
    PANEL_VALUE,"sample.wav",
    NULL);
  tsampling=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Sampling rate",
    Y,105,
    FIELD_TYPE,real_field,
    PANEL_REAL,"11025",
    NULL);
  tsamples=xv_create(fmain,TEXTFIELD,
    XV_LABEL,"Number of samples",
    X,255,
    Y,105,
    FIELD_TYPE,real_field,
    PANEL_REAL,"50000",
    NULL);
  bwrite=xv_create(fmain,BUTTON,
    XV_LABEL,"Write file",
    Y,135,
    NOTIFY_HANDLER,NotifyWrite,
    NULL);
  message1=xv_create(fmain,MESSAGE,
    X,99,
    Y,150,
    NULL);
  bquit=xv_create(fmain,BUTTON,
    XV_LABEL,"Quit",
    Y,150,
    NOTIFY_HANDLER,NotifyQuit,
    NULL);
  xv_main_loop(fmain);
  /* Exit */
  restorecrtmode();
}
