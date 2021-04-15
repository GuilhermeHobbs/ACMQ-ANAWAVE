/* Analysis and filtering of waveform files                                 */
/* By Antonio Carlos M. de Queiroz                                          */
/* Version 1.0 - 20/02/96                                                   */
/*
Waveform ~ok
Spectrogram ~ok.
Rectangles ok, Reading ok.
FFT ok.
FFT window ok.
Real FFT ok.
Play ok.
Record ok.
File save ok.
Working with djgpp.
Filtering ok.
Initialization ok. No more crashes due to log scales starting with 0.
21/12/96 Windowing.
25/12/96 Windowing with correct scaling.
1/1/97 Zoom in all modes in the FFT window. Correct file selection. Longer FFT.
2/1/97 Cursors.
8/1/97 Palette in color or b&w.
28/9/97 Corrected fft length limitation
14/4/99 Reads palette from file
23/02/2004 - Mixer disabled, recompiled with new mickey.

To do:
Allow shift step > fft length.
Better letters.
Partial replot in the FFT window if the panning keys are used.
Simpler reading routine.
Play section.
Phase in the FFT window?
Storage of FFT results.
Frequency in rd/s.
Vertical scale normalized to maximum=1.
Stereo files. Choose channel?
Cut and Paste
*/

#define debug(xxx) printf(xxx); getch(); close_window(active_w)

#define version "1.0c - 23/2/2004"

/* Undefine for complex FFT */
#define REALFFT

#ifdef REALFFT
#define GMAX 16384
#else
#define GMAX 8192
#endif
#define first_color 16
#define last_color 255
#define MINAMP 0.01
#define LOGMIN 0.00001
#define BMAX 30

#include <gppconio.h>
#include <libbcc.h>
#include <dpmi.h>
#define coreleft _go32_dpmi_remaining_virtual_memory

#include <string.h>
#include <math.h>
#include <ctype.h>
#include <dos.h>
#include <dir.h>
#include <process.h>
#include "xview.h"
#include "xv_var.h"
#include "extkeys.h"
#include "mickey.h"

/* Interface access */
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
#define get_string(xx) xx->v.ssetting.item_setting[xx->v.ssetting.sel_setting-1]
#define SECONDS (stime->v.ssetting.sel_setting==2)
#define LINEAR (slog->v.ssetting.sel_setting==1)
#define SQRT (slog->v.ssetting.sel_setting==2)
#define LOG (slog->v.ssetting.sel_setting==3)
#define COLOR_PALETTE (spalette->v.ssetting.sel_setting==1)
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
  l_int nSamplesPerSec; 	/* sampling frequency (Hz) */
  l_int nAvgBytesPerSec;	/* = nBlockAlign*nSamplesPerSec */
  s_int nBlockAlign;		/* = wBitsPerSample/8 * nChannels */
  s_int wBitsPerSample;		/* 16 or 8 */
  char data[4];		        /* "data" */
  l_int datasize;		/* size of the data block (bytes) */
} header;

unsigned s_int sample_size,channels,wave_start;
unsigned l_int wave_size;
FILE *wave;
s_int valid_file=0, valid_wave=0, valid_spectrum=0, valid_fft=0;
double ax,bx,ay,by,as,bs,af,bf,aa,ba,sampling_rate,avf,bvf,ahf,bhf;
s_int fftorder=7;
s_int grid=1;
double COS[GMAX], MSIN[GMAX], wmed;
s_int fxmin,fxmax,fymin,fymax,wxmin,wxmax,wymin,wymax,sxmin,sxmax,symin,symax;
char txt[256];

/* Declaration of the interface objects */
Xv_opaque menu2;
Xv_opaque fmessages,tty1;
Xv_opaque fio,tfilename,tmask,cdirectory,bread,bwrite,bplay;
Xv_opaque fpwave,tymax,tymin,txmin,txdelta,stime,
  bapplywave;
Xv_opaque fwave,cwave;
Xv_opaque fpspectrum,tfftsize,tshift,slog,swindow,tfmin,tfmax,sunit,
  tsmin,tsdelta,tamin,tamax,spalette,bapplyspectrum;
Xv_opaque fspectrum,cspectrum;
Xv_opaque fpfft,tfftsmin,tfftfmin,tfftfmax,tfftamin,tfftamax,bapplyfft,blistfft;
Xv_opaque ffft,cfft;
Xv_opaque ffilter,tfilter,tinfile,toutfile,tfxmin,tfxdelta,bapplyfilter;
Xv_opaque frecord,ssource,smode,ssampling,tsave,brecord,mrecord,ssize;
Xv_opaque fmixer,tmaster,tmic,tbalance,sagc,bapplymixer;

/* extern int _stklen=65000; */

/* To allow use of printf in the tty */
void tprintf(char* format,...)
{
  char teste[256];
  va_list paramlist;

  va_start(paramlist,format);
  vsprintf(teste,format,paramlist);
  ttysw_output(tty1,teste);
}
#define printf tprintf

void apply_window(double *A)
{
int m,i;
double k;

switch (get_sel(swindow)) {
  case 1: /* rectangular*/ return;
  case 2: /* triangular */
    m=get_int(tfftsize);
    for (i=0; i<m; i++) A[i]*=((1-fabs(2*(double)i/(double)m-1)))*2;
    break;
  case 3: /* parabolic */
    m=get_int(tfftsize);
    for (i=0; i<m; i++) {
      k=fabs(2*(double)i/(double)m-1);
      A[i]*=(1-(k*k))/wmed;
    }
    break;
  case 4: /* raised cosine */
    m=get_int(tfftsize);
    for (i=0; i<m; i++) A[i]*=(1-cos(2*M_PI/(double)m*(double)i));
  }
}

void PlotFFT(Xv_opaque obj);

void ComputeSinCos(void)
{
  s_int j,fftsize;
  double wt;

  fftsize=1<<fftorder;
  for (j=0; j<fftsize; j++) {
    wt=2*M_PI*j/fftsize;
    COS[j]=cos(wt);
    MSIN[j]=-sin(wt);
  }
}

s_int Lim(double x,double a,double b,s_int xmin,s_int xmax)
{
  double t;
  t=a*x+b;
  if (t>xmax) t=xmax;
  else if (t<xmin) t=xmin;
  return (s_int)t;
}

void close_dialogs(void)
{
  close_window(fpwave);
  close_window(fpspectrum);
  close_window(fpfft);
}

void DrawScales(double x1,double x2,double y1,double y2,
  s_int xmin, s_int xmax, s_int ymin, s_int ymax,
  s_int grid, s_int ylog)
{
  double ax,bx,ay,by,t1,t2;
  s_int i;
  char STR1[256];

  if (ylog==3) {
    ay=(ymax-ymin)/log(y1/y2);
    by=ymax-ay*log(y1);
  }
  else if (ylog==2) {
    ay=(ymax-ymin)/(sqrt(y1)-sqrt(y2));
    by=ymax-ay*sqrt(y1);
  }
  else {
    ay=(ymax-ymin)/(y1-y2);
    by=ymax-ay*y1;
  }
  ax=(xmax-xmin)/(x2-x1);
  bx=xmin-ax*x1;
  if (grid && x2>x1 && y2>y1) {
    setlinestyle(DOTTED_LINE,0,NORM_WIDTH);
    setcolor(LIGHTGRAY);
    if ((ylog==3)&&y2-y1>y1)
      t1=y1;
    else
      t1=y2-y1;
    t1=exp(log(10.0)*floor(log(t1)/log(10.0)-0.499999+0.5));
    t2=floor(y1/t1+1)*t1;
    while (t2<y2) {
      if (ylog!=3) {
        if (ylog==1) i=(s_int)floor(ay*t2+by+0.5);
        else i=(s_int)floor(ay*sqrt(t2)+by+0.5);
	line(xmin,i,xmax,i);
	t2+=t1;
	continue;
      }
      if ((s_int)floor(t2/t1+0.5)==10) {
	t1=10*t1;
	setcolor(YELLOW);
      }
      i=(s_int)floor(ay*log(t2)+by+0.5);
      line(xmin,i,xmax,i);
      t2+=t1;
      setcolor(LIGHTGRAY);
    }
    t1=x2-x1;
    t1=exp(log(10.0)*floor(log(t1)/log(10.0)));
    t2=floor(x1/t1+1)*t1;
    while (t2<x2) {
      i=(s_int)floor(ax*t2+bx+0.5);
      line(i,ymin,i,ymax);
      t2+=t1;
    }
  }
  setlinestyle(SOLID_LINE,0,NORM_WIDTH);
  setcolor(WHITE);
  line(0,ymax,xmax,ymax);
  line(xmin,ymin,xmin,getmaxy());
  sprintf(txt,"%8g",y2);
  outtextxy(0,ymin,txt);
  sprintf(txt,"%8g",y1);
  outtextxy(0,ymax-8,txt);
  sprintf(txt,"%8g",x1);
  while (txt[0]==' ') {
    strcpy(STR1,txt+1);
    strcpy(txt,STR1);
  }
  outtextxy(xmin+2,ymax+2,txt);
  sprintf(txt,"%8g",x2);
  outtextxy(xmax-64,ymax+2,txt);
  outtextxy((xmin+xmax-8*strlen(get_value(tfilename)))/2,ymax+2,get_value(tfilename));
}

double Imag;

double Cmult(double x1,double x2,double y1,double y2)
{
  Imag=x1*y2+x2*y1;
  return x1*y1-x2*y2;
}

#include "fht.c"

/* Callback procedures */

void InvalidateWave(Xv_opaque obj)
{
  valid_wave=0;
}

void InvalidateFFT(Xv_opaque obj)
{
  valid_fft=0;
  if ((SQRT || LOG)&& get_real(tfftamin)<MINAMP) set_real(tfftamin,MINAMP);
}

void InvalidateSpectrum(Xv_opaque obj)
{
  valid_spectrum=0;
  if ((SQRT || LOG)&& get_real(tamin)<MINAMP) {
    set_real(tamin,MINAMP);
    set_real(tfftamin,MINAMP);
  }
}

void ReadFile(Xv_opaque obj)
{
  /* Notify handler for bread */
  /* .WAV header reader adapted from code by David Welch, 1995 */
  char gstring[80];
  unsigned short ra;
  unsigned short rb;
  unsigned short rc;
  unsigned long la;
  unsigned short fmttag;
  unsigned long ckid;
  unsigned long cksize;
  unsigned long filesize;
  unsigned long nextseek;
  s_int k;

  fmttag=0;
  if((wave=fopen(get_value(tfilename),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tfilename));
    return;
  }
  printf("Reading file %s:\n",get_value(tfilename));
  fread(&la,1,4,wave);
  if(la!=0x46464952) {
    printf("Not a RIFF format file\n");
    goto End;
  }
  close_window(fio);
  fread(&filesize,1,4,wave);
  filesize+=ftell(wave);
  printf("File size: %lu bytes\n",filesize);
  fread(&la,1,4,wave);
  if(la!=0x45564157)
  {
    printf("Not a WAVE file\n");
    goto End;
  }
  nextseek=ftell(wave);
  while(nextseek<filesize) {
    fseek(wave,nextseek,0);
    fread(&ckid,1,4,wave);
    fread(&cksize,1,4,wave);
    nextseek=cksize+ftell(wave);
    memcpy(gstring,(void *)&ckid,4); gstring[4]=0;
    printf("Chunk ID: %s; Size %lu\n",gstring,cksize);
    switch(ckid) {
      case 0x20746D66: //Format Chunk
        fread(&fmttag,1,2,wave);
        printf("  Format tag %04Xh:  ",fmttag);
        switch(fmttag) {
          default:
          case 0x0000: printf("WAVE_FORMAT_UNKNOWN\n"); break;
          case 0x0001: printf("WAVE_FORMAT_PCM\n"); break;
          case 0x0002: printf("WAVE_FORMAT_ADPCM\n"); break;
          case 0x0005: printf("WAVE_FORMAT_IBM_CVSD\n"); break;
          case 0x0006: printf("WAVE_FORMAT_ALAW\n"); break;
          case 0x0007: printf("WAVE_FORMAT_MULAW\n"); break;
          case 0x0010: printf("WAVE_FORMAT_OKI_ADPCM\n"); break;
          case 0x0011: printf("WAVE_FORMAT_DVI_ADPCM or WAVE_FORMAT_IMA_ADPCM\n"); break;
          case 0x0015: printf("WAVE_FORMAT_DIGISTD\n"); break;
          case 0x0016: printf("WAVE_FORMAT_DIGIFIX\n"); break;
          case 0x0020: printf("WAVE_FORMAT_YAMAHA_ADPCM\n"); break;
          case 0x0021: printf("WAVE_FORMAT_SONARC\n"); break;
          case 0x0022: printf("WAVE_FORMAT_DSPGROUP_TRUESPEECH\n"); break;
          case 0x0023: printf("WAVE_FORMAT_ECHOSC1\n"); break;
          case 0x0024: printf("WAVE_FORMAT_AUDIOFILE_AF18\n"); break;
          case 0x0101: printf("IBM_FORMAT_MULAW\n"); break;
          case 0x0102: printf("IBM_FORMAT_ALAW\n"); break;
          case 0x0103: printf("IBM_FORMAT_ADPCM\n"); break;
          case 0x0200: printf("WAVE_FORMAT_CREATIVE_ADPCM\n"); break;
        }
        fread(&ra,1,2,wave);
        channels=ra;
        printf("  %u channels\n",ra);
        fread(&la,1,4,wave);
        sampling_rate=(double)la;
        printf("  %lu samples per second\n",la);
        fread(&la,1,4,wave);
        printf("  %lu average bytes per second\n",la);
        fread(&ra,1,2,wave);
        printf("  Block align %u\n",ra);
        fread(&ra,1,2,wave);
        sample_size=ra/8;
        printf("  %u bits per sample\n",ra);
        if(fmttag!=0x0001) {
          fread(&ra,1,2,wave);
          printf("  %u bytes of extra information\n");
          switch(fmttag) {
            case 0x0002:    //MS ADPCM?
              fread(&ra,1,2,wave);
              printf("  %u samples per block\n",ra);
              fread(&ra,1,2,wave);
              printf("  %u coefficients:\n",ra);
              for(la=0;la<ra;la++) {
                fread(&rb,1,2,wave);
                fread(&rc,1,2,wave);
                printf("    %u  %u\n",rb,rc);
              }
              break;
            case 0x0011: //WAVE_FORMAT_DVI_ADPCM
            case 0x0022: //WAVE_FORMAT_DSPGROUP_TRUESPEECH
              fread(&ra,1,2,wave);
              printf("  %u samples per block\n",ra);
              break;
            case 0x0021:  // WAVE_FORMAT_SONARC
              fread(&ra,1,2,wave);
              printf("  %u CompType\n",ra);
              break;
            case 0x0200: //WAVE_FORMAT_CREATIVE_ADPCM
              fread(&ra,1,2,wave);
              printf("  Revision %u\n",ra);
              break;
          }
        }
        break;
      case 0x61746164: //Data Chunk
        if (!fmttag) {
          printf("Format chunk not defined before data chunk\n");
          goto End;
        }
        wave_start=ftell(wave);
        wave_size=cksize/sample_size;
        break;
      default:
        printf("Warning: Unknown chunk\n");
    }
  }
  if (fmttag!=0x0001) {
    printf("Unsupported format\n");
    goto End;
  }
  if (channels!=1) {
    printf("Stereo files are not supported\n");
    goto End;
  }
  set_real(txmin,1);
  set_real(txdelta,wave_size-1);
  set_real(tsmin,1);
  set_real(tsdelta,wave_size-1);
  set_real(tfxmin,1);
  set_real(tfxdelta,wave_size-1);
  if (sample_size==1) k=127; else k=32767;
  set_real(tymin,-k-1);
  set_real(tymax,k);
  set_real(tamin,0);
  set_real(tamax,k);
  set_real(tfftamin,0);
  set_real(tfftamax,k);
  valid_file=1;
  close_window(ffilter);
  printf("The file %s is valid\n",get_value(tfilename));
  InvalidateSpectrum(NULL);
  InvalidateWave(NULL);
  InvalidateFFT(NULL);
  open_window(fpwave);
  open_window(fpspectrum);
  /*open_window(fpfft);*/
 End:
  fclose(wave);
}

void WriteFile(Xv_opaque obj)
{
  /* Notify handler for bwrite */
  open_window(frecord);
}

void NotifyWave(Xv_opaque obj)
{
  /* Notify handler for cwave */
  unsigned l_int s,s1,s2;
  s_int sample,xx,yy,xxa,c0,c1,c2,maxyy,minyy;
  double x1,x2;

  if (!valid_file) return;
  if (get_real(txmin)<1) set_real(txmin,1);
  if (get_real(txdelta)<1) set_real(txdelta,1);
  if (get_real(txdelta)>wave_size-1) set_real(txdelta,wave_size-1);
  if (get_real(txmin)+get_real(txdelta)>wave_size) set_real(txmin,wave_size-get_real(txdelta));
  wxmin=65; wymin=3; wxmax=get_dx(cwave)-5; wymax=get_dy(cwave)-11;
  ay=(wymax-wymin)/(get_real(tymin)-get_real(tymax));
  by=wymax-ay*get_real(tymin);
  ax=(double)(wxmax-wxmin)/get_real(txdelta);
  bx=wxmin-ax*get_real(txmin);
  setfillstyle(SOLID_FILL,BLACK);
  bar(0,0,get_dx(cwave)-2,get_dy(cwave)-2);
  x1=get_real(txmin); x2=x1+get_real(txdelta);
  if (SECONDS) {x1=(x1-1)/sampling_rate; x2=(x2-1)/sampling_rate;}
  DrawScales(x1,x2,get_real(tymin),get_real(tymax),wxmin,wxmax,wymin,wymax,grid,1);
  if ((wave=fopen(get_value(tfilename),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tfilename));
    return;
  }
  s1=(unsigned l_int)get_real(txmin);
  s2=s1+(unsigned l_int)get_real(txdelta);
  fseek(wave,wave_start+(s1-1)*sample_size,0);
  setcolor(WHITE);
  xxa=-10; maxyy=-10; minyy=10000;
  for (s=s1; s<=s2; s++) {
    fread(&sample,sample_size,1,wave);
    if (sample_size==1) sample=(sample & 255)-128;
    xx=Lim((double)s,ax,bx,wxmin,wxmax);
    yy=Lim((double)sample,ay,by,wymin,wymax);
    if (s==s1) moveto(xx,yy);
    else {
      if ((c0=xx)!=xxa) {minyy=maxyy=yy; xxa=xx;}
      else {
	if ((c1=yy)>maxyy) maxyy=yy;
	if ((c2=yy)<minyy) minyy=yy;
      }
      if (c0 || c1 || c2) {
	lineto(xx,yy);
	if (mkbhit()) goto Fim;
      }
    }
  }
 Fim:
  fclose(wave);
  valid_wave=1;
}

void EventWave(Xv_opaque obj)
{
  /* Event handler for cwave */
  static s_int xi=0, yi=0, xf=0, yf=0;
  static s_int selecao=0;
  double t;

  if (!valid_wave) return;
  if (selecao) {
    setwritemode(XOR_PUT);
    setcolor(WHITE);
    if (ie_code==LOC_DRAG) {
      rectangle(xi,yi,xf,yf);
      xf=ie_locx-1;
      yf=ie_locy-1;
      rectangle(xi,yi,xf,yf);
      return;
    }
    if (ie_shiftcode==0) {
      selecao=0;
      rectangle(xi,yi,xf,yf);
      if (xi>=xf||yi>=yf) return;
      t=(xi-bx)/ax;
      if (t<1) t=1;
      set_real(txmin,(unsigned l_int)t);
      set_real(txdelta,(unsigned l_int)((xf-xi)/ax));
      set_real(tymin,(yf-by)/ay);
      set_real(tymax,(yi-by)/ay);
    }
    setwritemode(COPY_PUT);
    goto Recalcular;
  }

  if (ie_code<256) ie_code=toupper(ie_code);
  switch (ie_code) {

    case MS_MIDDLE:
    case 'S':
      t=(ie_locx-1-bx)/ax;
      if (t<0) t=0;
      set_real(tfftsmin,(unsigned l_int)t);
      PlotFFT(NULL);
      return;

    case 'R':
      set_real(txdelta,get_real(txdelta)*2);
      goto Recalcular;

    case 'A':
      set_real(txdelta,floor(get_real(txdelta)/2));
      goto Recalcular;

    case '>':
    case '.':
      txmin->v.stextfield.panel_real+=(unsigned l_int)(get_real(txdelta)/4);
      goto Recalcular;

    case '<':
    case ',':
      txmin->v.stextfield.panel_real-=(unsigned l_int)(get_real(txdelta)/4);
      goto Recalcular;

    case 'G':
      grid=!grid;
      goto Recalcular;

    case KUPARROW:
      t=(get_real(tymax)-get_real(tymin))/4;
      tymin->v.stextfield.panel_real+=t;
      tymax->v.stextfield.panel_real+=t;
      goto Recalcular;

    case KDOWNARROW:
      t=(get_real(tymax)-get_real(tymin))/4;
      tymin->v.stextfield.panel_real-=t;
      tymax->v.stextfield.panel_real-=t;
      goto Recalcular;

    case '-':
      set_real(tymax,2*get_real(tymax)-get_real(tymin));
      goto Recalcular;

    case '+':
      set_real(tymax,(get_real(tymax)+get_real(tymin))/2);
      goto Recalcular;

    case MS_LEFT:
      xi=ie_locx-1;
      yi=ie_locy-1;
      xf=xi;
      yf=yi;
      selecao=1;
      setcolor(WHITE);
      rectangle(xi,yi,xf,yf);
      return;

    case ' ':
      t=(ie_locx-1-bx)/ax;
      if (SECONDS) t=(t-1)/sampling_rate;
      sprintf(txt,"H=%g V=%g",t,(ie_locy-1-by)/ay);
      setfillstyle(SOLID_FILL,BLACK);
      bar(wxmin+1,wymin-1,wxmin+240,wymin+8);
      outtextxy(wxmin+2,wymin,txt);

    default: return;
  }
Recalcular:
  NotifyWave(NULL);
}

void PlotWave(Xv_opaque obj)
{
  /* Notify handler for bapplywave */
  s_int k;

  close_dialogs();
  k=fwave->v.sframe.mapped;
  open_window(fwave);
  if (k && xv_ok) NotifyWave(NULL);
}

void ChangeWindow(Xv_opaque obj)
{
  int i,m;
  double k;
  /* Apenas o valor para a janela parabolica e usado. Os demais sao fixos */
  m=get_int(tfftsize); wmed=0;
  switch (get_sel(swindow)) {
    case 1: /* rectangular*/ wmed=m; break;
    case 2: /* triangular */
      for (i=0; i<m; i++) wmed+=(1-fabs(2*(double)i/(double)m-1));
    break;
    case 3: /* parabolic */
      for (i=0; i<m; i++) {
        k=fabs(2*(double)i/(double)m-1);
        wmed+=(1-(k*k));
      }
      break;
    case 4: /* raised cosine */
      for (i=0; i<m; i++) wmed+=0.5*(1-cos(2*M_PI/(double)m*(double)i));
  }
  wmed/=m;
  /*
  printf("Window medium value = %f10.8\n",wmed);
  */
  InvalidateSpectrum(NULL);
}

void NFFTLength(Xv_opaque obj)
{
  /* Notify handler for tfftsize */
  l_int k;

  k=1;
  while (k<get_int(tfftsize)) k*=2;
  if (k>GMAX) k=GMAX;
  set_int(tfftsize,k);
  update(tfftsize);
  fftorder=0; k=1;
  while(k<get_int(tfftsize)) {
    fftorder++;
    k*=2;
  }
  ComputeSinCos();
  ChangeWindow(NULL);
}

void NotifySpectrum(Xv_opaque obj)
{
  /* notify handler for cspectrum */
  unsigned l_int s,s1,s2,sm;
  s_int sample,yy,yya,j,k,ds,sb,rest,n,T[GMAX];
  double df,A[GMAX],t,max,x1,x2;
#ifndef REALFFT
  double B[GMAX];
#endif
  if (!valid_file) return;
  if (get_int(tshift)>get_int(tfftsize)) set_int(tshift,get_int(tfftsize));
  if (get_real(tsmin)<1) set_real(tsmin,1);
  if (get_real(tsdelta)<1) set_real(tsdelta,1);
  if (get_real(tsdelta)>wave_size-1) set_real(tsdelta,wave_size-1);
  if (get_real(tsmin)+get_real(tsdelta)>wave_size) set_real(tsmin,wave_size-get_real(tsdelta));
  sxmin=65; symin=3; sxmax=get_dx(cspectrum)-5; symax=get_dy(cspectrum)-11;
  af=(symax-symin)/(get_real(tfmin)-get_real(tfmax));
  bf=symax-af*get_real(tfmin);
  as=(double)(sxmax-sxmin)/get_real(tsdelta);
  bs=sxmin-as*get_real(tsmin);
  x1=get_real(tamin);
  x2=get_real(tamax);
  if (SQRT) {
    aa=(last_color-first_color)/(sqrt(x2)-sqrt(x1));
    ba=first_color-aa*sqrt(x1);
  }
  else if (LOG) {
    aa=(last_color-first_color)/(log(x2)-log(x1));
    ba=first_color-aa*log(x1);
  }
  else {
    aa=(last_color-first_color)/(x2-x1);
    ba=first_color-aa*x1;
  }
  setfillstyle(SOLID_FILL,BLACK);
  bar(0,0,get_dx(cspectrum)-2,get_dy(cspectrum)-2);
  if ((wave=fopen(get_value(tfilename),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tfilename));
    return;
  }
  s1=(unsigned l_int)get_real(tsmin);
  s2=s1+(unsigned l_int)get_real(tsdelta);
  fseek(wave,wave_start+(s1-1)*sample_size,0);
  df=af*sampling_rate/(double)(get_int(tfftsize));
  ds=get_int(tshift);
  rest=get_int(tfftsize)-ds;
  max=0;
  x1=get_real(tsmin); x2=x1+get_real(tsdelta);
  if (SECONDS) {x1=(x1-1)/sampling_rate; x2=(x2-1)/sampling_rate;}
  DrawScales(x1,x2,get_real(tfmin),get_real(tfmax),sxmin,sxmax,symin,symax,grid,1);
  n=get_int(tfftsize); k=n/2;
  for (s=s1; s<=s2; s+=ds) {
    for (sb=0; sb<n; sb++) {
      if ((s==s1) || sb>=rest) {
	if (fread(&sample,sample_size,1,wave)==1) {
	  if (sample_size==1) sample=(sample & 255)-128;
	}
	else sample=0;
      }
      else sample=T[sb+ds];
      T[sb]=sample;
      A[sb]=(double)sample;
#ifndef REALFFT
      B[sb]=0.0;
#endif
    }
    apply_window(A);
#ifdef REALFFT
    realfft(n,A);
    A[0]=fabs(A[0])/n;
    for (j=1; j<k; j++) A[j]=sqrt(A[j]*A[j]+A[n-j]*A[n-j])/k;
#else
    fft(n,A,B);
    A[0]=fabs(A[0])/n;
    for (j=1; j<k; j++) A[j]=sqrt(A[j]*A[j]+B[j]*B[j])/k;
#endif
    t=-1000;
    yya=0;
    for (j=0; j<k; j++) {
      if (A[j]>max) {max=A[j]; sm=s;}
      if (SQRT) A[j]=sqrt(A[j]);
      else if (LOG) A[j]=log(A[j]+LOGMIN);
      yy=Lim((double)j,df,bf,symin,symax);
      if ((yy!=yya) || (A[j]>t)) {
        setfillstyle(SOLID_FILL,Lim(A[j],aa,ba,first_color,last_color));
        bar(Lim((double)s,as,bs,sxmin,sxmax),
          yy,
          Lim((double)(s+ds),as,bs,sxmin,sxmax),
          Lim((double)(j+1),df,bf,symin,symax));
        yya=yy;
        t=A[j];
      }
    }
    if (mkbhit()) goto Fim;
  }
 Fim:
  DrawScales(x1,x2,get_real(tfmin),get_real(tfmax),sxmin,sxmax,symin,symax,grid,1);
  fclose(wave);
  /*
  if (fmessages->v.sframe.mapped) {
    printf("Spectrogram:\n  Maximum amplitude: %g, at sample %lu\n",max,sm);
  }
  */
  valid_spectrum=1;
}

void EventSpectrum(Xv_opaque obj)
{
  /* Event handler for cspectrum */
  static s_int xi=0, yi=0, xf=0, yf=0;
  static s_int selecao=0;
  double t;

  if (!valid_spectrum) return;
  if (selecao) {
    setwritemode(XOR_PUT);
    setcolor(WHITE);
    if (ie_code==LOC_DRAG) {
      rectangle(xi,yi,xf,yf);
      xf=ie_locx-1;
      yf=ie_locy-1;
      rectangle(xi,yi,xf,yf);
      return;
    }
    if (ie_shiftcode==0) {
      selecao=0;
      rectangle(xi,yi,xf,yf);
      if (xi>=xf||yi>=yf) return;
      t=(xi-bs)/as;
      if (t<1) t=1;
      set_real(tsmin,(unsigned l_int)t);
      set_real(tsdelta,(unsigned l_int)((xf-xi)/as));
      set_real(tfmin,(yf-bf)/af);
      set_real(tfmax,(yi-bf)/af);
    }
    setwritemode(COPY_PUT);
    goto Recalcular;
  }
  if (ie_code<256) ie_code=toupper(ie_code);
  switch (ie_code) {

    case MS_MIDDLE:
    case 'S':
      t=(ie_locx-1-bs)/as;
      if (t<0) t=0;
      set_real(tfftsmin,(unsigned l_int)t);
      set_real(tfftamin,get_real(tamin));
      set_real(tfftamax,get_real(tamax));
      set_real(tfftfmin,get_real(tfmin));
      set_real(tfftfmax,get_real(tfmax));
      PlotFFT(NULL);
      return;

    case 'R':
      set_real(tsdelta,get_real(tsdelta)*2);
      goto Recalcular;

    case 'A':
      set_real(tsdelta,floor(get_real(tsdelta)/2));
      goto Recalcular;

    case '>':
    case '.':
      tsmin->v.stextfield.panel_real+=(unsigned l_int)(get_real(tsdelta)/4);
      goto Recalcular;

    case '<':
    case ',':
      tsmin->v.stextfield.panel_real-=(unsigned l_int)(get_real(tsdelta)/4);
      goto Recalcular;

    case 'G':
      grid=!grid;
      goto Recalcular;

    case KUPARROW:
      t=(get_real(tfmax)-get_real(tfmin))/4;
      tfmin->v.stextfield.panel_real+=t;
      tfmax->v.stextfield.panel_real+=t;
      goto Recalcular;

    case KDOWNARROW:
      t=(get_real(tfmax)-get_real(tfmin))/4;
      tfmin->v.stextfield.panel_real-=t;
      tfmax->v.stextfield.panel_real-=t;
      goto Recalcular;

    case '-':
      set_real(tfmax,2*get_real(tfmax)-get_real(tfmin));
      goto Recalcular;

    case '+':
      set_real(tfmax,(get_real(tfmax)+get_real(tfmin))/2);
      goto Recalcular;

    case MS_LEFT:
      xi=ie_locx-1;
      yi=ie_locy-1;
      xf=xi;
      yf=yi;
      selecao=1;
      setcolor(WHITE);
      rectangle(xi,yi,xf,yf);
      return;

    case ' ':
      t=(ie_locx-1-bs)/as;
      if (SECONDS) t=(t-1)/sampling_rate;
      sprintf(txt,"H=%g F=%g Hz",t,(ie_locy-1-bf)/af);
      setfillstyle(SOLID_FILL,BLACK);
      bar(sxmin+1,symin-1,sxmin+240,symin+8);
      outtextxy(sxmin+1,symin,txt);

    default: return;
  }
Recalcular:
  NotifySpectrum(NULL);
}

void PlotSpectrum(Xv_opaque obj)
{
  /* Notify handler for bapplyspectrum */
  s_int k;

  close_dialogs();
  k=fspectrum->v.sframe.mapped;
  open_window(fspectrum);
  if (k && xv_ok) NotifySpectrum(NULL);
}

void NotifyFFT(Xv_opaque obj)
{
  /* notify handler for cfft */
  unsigned l_int s1;
  s_int sample,xx,yy,j,k,sb,fm,n;
  double A[GMAX],max,t;
#ifndef REALFFT
  double B[GMAX];
#endif
  char tmp[30];

  if (!valid_file) return;
  if (get_real(tfftsmin)<1) set_real(tfftsmin,1);
  if (get_real(tfftsmin)>wave_size) set_real(tfftsmin,wave_size-1);
  fxmin=65; fymin=3; fxmax=get_dx(cfft)-5; fymax=get_dy(cfft)-11;
  ahf=(fxmax-fxmin)/(get_real(tfftfmax)-get_real(tfftfmin));
  bhf=fxmax-ahf*get_real(tfftfmax);
  if (SQRT) {
    avf=(fymin-fymax)/(sqrt(get_real(tfftamax))-sqrt(get_real(tfftamin)));
    bvf=fymax-avf*sqrt(get_real(tfftamin));
  }
  else if (LOG) {
    avf=(fymin-fymax)/(log(get_real(tfftamax))-log(get_real(tfftamin)));
    bvf=fymax-avf*log(get_real(tfftamin));
  }
  else {
    avf=(fymin-fymax)/(get_real(tfftamax)-get_real(tfftamin));
    bvf=fymax-avf*get_real(tfftamin);
  }
  setfillstyle(SOLID_FILL,BLACK);
  bar(0,0,get_dx(cfft)-2,get_dy(cfft)-2);
  if ((wave=fopen(get_value(tfilename),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tfilename));
    return;
  }
  t=(double)(fymax-fymin)/(double)(first_color-last_color);
  for (j=first_color; j<=last_color; j++) {
    setcolor(j);
    yy=(s_int)(t*(j-last_color)+fymin);
    line(3,yy,20,yy);
    line(3,yy-1,20,yy-1);
  }
  DrawScales(get_real(tfftfmin),get_real(tfftfmax),get_real(tfftamin),get_real(tfftamax),fxmin,fxmax,fymin,fymax,grid,get_sel(slog));
  s1=(unsigned l_int)get_real(tfftsmin);
  sprintf(tmp,"Spectrum at samples %lu-%lu",s1,s1+(unsigned l_int)(get_int(tfftsize)-1));
  xv_set(ffft,tmp);
  setcolor(WHITE);
  n=get_int(tfftsize);
  k=n/2;
  fseek(wave,wave_start+(s1-1)*sample_size,0);
  for (sb=0; sb<n; sb++) {
    if (fread(&sample,sample_size,1,wave)==1) {
      if (sample_size==1) sample=(sample & 255)-128;
      A[sb]=(double)sample;
    }
    else A[sb]=0.0;
#ifndef REALFFT
    B[sb]=0.0;
#endif
  }
  apply_window(A);
#ifdef REALFFT
  realfft(n,A);
  A[0]=fabs(A[0])/n;
  for (j=1; j<k; j++) A[j]=sqrt(A[j]*A[j]+A[n-j]*A[n-j])/k;
#else
  fft(n,A,B);
  A[0]=fabs(A[0])/n;
  for (j=1; j<k; j++) A[j]=sqrt(A[j]*A[j]+B[j]*B[j])/k;
#endif
  max=-1;
  for (j=0; j<k; j++) {
    if (A[j]>max) {max=A[j]; fm=j;}
    if (SQRT) A[j]=sqrt(A[j]);
    else if (LOG) A[j]=log(A[j]+LOGMIN);
    xx=Lim((double)j*sampling_rate/get_int(tfftsize),ahf,bhf,fxmin,fxmax);
    yy=Lim(A[j],avf,bvf,fymin,fymax);
    if (j==0) moveto(xx,yy);
    else lineto(xx,yy);
    if (mkbhit()) goto Fim;
  }
 Fim:
  fclose(wave);
  /*
  if (fmessages->v.sframe.mapped) {
    printf("Spectrum of samples %lu - %lu\n  Max. amplitude: %g, at %g Hz\n",s1,s1+(unsigned l_int)(get_int(tfftsize)-1),max,fm*sampling_rate/get_int(tfftsize));
  }
  */
  valid_fft=1;
}

void EventFFT(Xv_opaque obj)
{
  /* Event handler for cfft */

  static s_int xi=0, yi=0, xf=0, yf=0;
  static s_int selecao=0;
  double t1,t2;

  if (!valid_fft) return;
  if (selecao) {
    setwritemode(XOR_PUT);
    setcolor(WHITE);
    if (ie_code==LOC_DRAG) {
      rectangle(xi,yi,xf,yf);
      xf=ie_locx-1;
      yf=ie_locy-1;
      rectangle(xi,yi,xf,yf);
      return;
    }
    if (ie_shiftcode==0) {
      selecao=0;
      rectangle(xi,yi,xf,yf);
      if (xi>=xf||yi>=yf) return;
      set_real(tfftfmax,(xf-bhf)/ahf);
      set_real(tfftfmin,(xi-bhf)/ahf);
      t1=(yi-bvf)/avf;
      t2=(yf-bvf)/avf;
      if (LINEAR) {
        set_real(tfftamax,t1);
        set_real(tfftamin,t2);
      }
      else if (SQRT) {
        set_real(tfftamax,t1*t1);
        set_real(tfftamin,t2*t2);
      }
      else {
        set_real(tfftamax,exp(t1)-LOGMIN);
        set_real(tfftamin,exp(t2)-LOGMIN);
      }
    }
    setwritemode(COPY_PUT);
    goto Recalcular;
  }
  if (ie_code<256) ie_code=toupper(ie_code);
  switch (ie_code) {

    case 'R':
      set_real(tfftfmax,get_real(tfftfmin)+(get_real(tfftfmax)-get_real(tfftfmin))*2);
      goto Recalcular;

    case 'A':
      set_real(tfftfmax,get_real(tfftfmin)+(get_real(tfftfmax)-get_real(tfftfmin))/2);
      goto Recalcular;

    case '>':
    case '.':
      t1=(get_real(tfftfmax)-get_real(tfftfmin))/4;
      tfftfmin->v.stextfield.panel_real+=t1;
      tfftfmax->v.stextfield.panel_real+=t1;
      goto Recalcular;

    case '<':
    case ',':
      t1=(get_real(tfftfmax)-get_real(tfftfmin))/4;
      tfftfmin->v.stextfield.panel_real-=t1;
      tfftfmax->v.stextfield.panel_real-=t1;
      goto Recalcular;

    case '}':
      tfftsmin->v.stextfield.panel_real+=1;
      goto Recalcular;

    case '{':
      tfftsmin->v.stextfield.panel_real-=1;
      goto Recalcular;

    case ']':
      tfftsmin->v.stextfield.panel_real+=get_int(tshift);
      goto Recalcular;

    case '[':
      tfftsmin->v.stextfield.panel_real-=get_int(tshift);
      goto Recalcular;

    case 'G':
      grid=!grid;
      goto Recalcular;

    case KUPARROW:
      if (LOG) {
	t1=sqrt(sqrt(get_real(tfftamax)/get_real(tfftamin)));
	tfftamin->v.stextfield.panel_real*=t1;
	tfftamax->v.stextfield.panel_real*=t1;
      }
      else if (SQRT) {
	t1=(sqrt(get_real(tfftamax))-sqrt(get_real(tfftamin)))/4;
	t2=sqrt(get_real(tfftamin))+t1;
	set_real(tfftamin,t2*t2);
	t2=sqrt(get_real(tfftamax))+t1;
	set_real(tfftamax,t2*t2);
      }
      else {
	t1=(get_real(tfftamax)-get_real(tfftamin))/4;
	tfftamin->v.stextfield.panel_real+=t1;
	tfftamax->v.stextfield.panel_real+=t1;
      }
      goto Recalcular;

    case KDOWNARROW:
      if (LOG) {
	t1=sqrt(sqrt(get_real(tfftamax)/get_real(tfftamin)));
	tfftamin->v.stextfield.panel_real/=t1;
	tfftamax->v.stextfield.panel_real/=t1;
      }
      else if (SQRT) {
	t2=sqrt(get_real(tfftamin));
	t1=(sqrt(get_real(tfftamax))-t2)/4;
	if (t1<t2) {
	  t2=t2-t1;
	  set_real(tfftamin,t2*t2);
	  t2=sqrt(get_real(tfftamax))-t1;
	}
	else {
	  set_real(tfftamin,MINAMP);
	  t2=sqrt(MINAMP)+4*t1;
	}
      set_real(tfftamax,t2*t2);
      }
      else {
	t1=(get_real(tfftamax)-get_real(tfftamin))/4;
	tfftamin->v.stextfield.panel_real-=t1;
	tfftamax->v.stextfield.panel_real-=t1;
      }
      goto Recalcular;

    case '-':
      if (LOG) {
	t1=get_real(tfftamax)/get_real(tfftamin);
	set_real(tfftamax,get_real(tfftamin)*t1*t1);
      }
      else if (SQRT) {
	t1=2*sqrt(get_real(tfftamax))-sqrt(get_real(tfftamin));
	set_real(tfftamax,t1*t1);
      }
      else
	set_real(tfftamax,2*get_real(tfftamax)-get_real(tfftamin));
      goto Recalcular;

    case '+':
      if (LOG)
	set_real(tfftamax,get_real(tfftamin)*sqrt(get_real(tfftamax)/get_real(tfftamin)));
      else if (SQRT) {
	t1=0.5*(sqrt(get_real(tfftamax))+sqrt(get_real(tfftamin)));
	set_real(tfftamax,t1*t1);
      }
      else
	set_real(tfftamax,(get_real(tfftamax)+get_real(tfftamin))/2);
      goto Recalcular;

    case MS_LEFT:
      xi=ie_locx-1;
      yi=ie_locy-1;
      xf=xi;
      yf=yi;
      selecao=1;
      setcolor(WHITE);
      rectangle(xi,yi,xf,yf);
      return;

    case ' ':
      t1=(ie_locx-1-bhf)/ahf;
      t2=(ie_locy-1-bvf)/avf;
      if (SQRT) t2*=t2;
      else if (LOG) t2=exp(t2)-LOGMIN;
      sprintf(txt,"F=%g Hz M=%g",t1,t2);
      setfillstyle(SOLID_FILL,BLACK);
      bar(fxmin+1,fymin-1,fxmin+240,fymin+8);
      outtextxy(fxmin+1,fymin,txt);

    default: return;
  }
Recalcular:
  NotifyFFT(NULL);
}

void PlotFFT(Xv_opaque obj)
{
  /* Notify handler for bapplyfft */
  s_int k;
  close_dialogs();
  k=ffft->v.sframe.mapped;
  open_window(ffft);
  if (k && xv_ok) NotifyFFT(NULL);
}

void ListFFT(Xv_opaque obj)
{
  /* notify handler for blistfft */
  unsigned l_int s1;
  s_int sample,j,sb,n;
  double A[GMAX];
#ifndef REALFFT
  double B[GMAX];
#endif
  if (!valid_file) return;
  if (get_real(tfftsmin)<1) set_real(tfftsmin,1);
  if (get_real(tfftsmin)>wave_size) set_real(tfftsmin,wave_size-1);
  if ((wave=fopen(get_value(tfilename),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tfilename));
    return;
  }
  s1=(unsigned l_int)get_real(tfftsmin);
  n=get_int(tfftsize);
  printf("Signal at samples %lu-%lu\n",s1,s1+(unsigned l_int)(n-1));
  fseek(wave,wave_start+(s1-1)*sample_size,0);
  for (sb=0; sb<n; sb++) {
    if (fread(&sample,sample_size,1,wave)==1) {
      if (sample_size==1) sample=(sample & 255)-128;
    }
    else sample=0;
    A[sb]=(double)sample;
#ifndef REALFFT
    B[sb]=0.0;
#endif
    printf("f[%d]=%d\n",sb,sample);
  }
  apply_window(A);
#ifdef REALFFT
  realfft(n,A);
  printf("Spectrum at samples %lu-%lu\n",s1,s1+(unsigned l_int)(n-1));
  for (j=0; j<n/2; j++)
    printf("F[%d]: %g %gj |F|=%g\n",j,
      A[j]/n,
      A[n-j]/n,
      sqrt(A[j]*A[j]+A[n-j]*A[n-j])/n);
#else
  fft(n,A,B);
  printf("Spectrum at samples %lu-%lu\n",s1,s1+(unsigned l_int)(n-1));
  for (j=0; j<n/2; j++)
    printf("F[%d]: %g %gj |F|=%g\n",j,
      A[j]/n,
      B[j]/n,
      sqrt(A[j]*A[j]+B[j]*B[j])/n);
#endif
  fclose(wave);
}

void PlayFile(Xv_opaque obj)
{
  char *blaster, path[100];
  s_int result;

  if ((blaster=getenv("SOUND"))==NULL) {
    printf("SOUND environment variable not found\n");
    return;
  }
  printf("Playing the file %s\nTouch ESC to interrupt\n",get_value(tfilename));
  sprintf(path,"%s/%s",blaster,"play");
  result=spawnl(P_WAIT,path,"1",get_value(tfilename),"/Q",NULL);
  if (result) printf("Some error occurred\n");
  printf("Free memory: %lu\n",coreleft());
}

void NotifyMixer(Xv_opaque obj)
{
  printf("Not yet implemented\n");
/*
  Some text is printed at the upper left corner if AGC is used,
  and the mouse moves to the lower right corner

  char *blaster, path[100], t2[100];
  s_int result;

  if ((blaster=getenv("SOUND"))==NULL) {
    printf("SOUND environment variable not found\n");
    return;
  }
  sprintf(path,"%s/%s",blaster,"sb16set");
  sprintf(t2,"/Q /MA:%d;%d /MIC:%d /AGC:%c ",
    get_int(tmaster),
    get_int(tbalance),
    get_int(tmic),
    get_sel(sagc)==1?'+':'-'
    );

  printf("Command: %s %s\n",path,t2);
  result=spawnl(P_WAIT,path,"1",t2,NULL);
  if (result) printf("Some error occurred\n");
  printf("Free memory: %lu\n",coreleft());
*/
}

void ProcessMenu2(Xv_opaque obj)
{
  /* Notify handler for menu2 */
  switch (obj->v.smenu.sel_menu) {
    case 1: /* Read/write signal */
      open_window(fio);
      break;
    case 2: /* Plot waveform */
      open_window(fpwave);
      break;
    case 3: /* Plot spectrogram */
      open_window(fpspectrum);
      break;
    case 4: /* Plot spectrum */
      open_window(fpfft);
      break;
    case 5: /* Filter signal */
      set_value(tinfile,get_value(tfilename));
      update(tinfile);
      break;
    case 6: /* Play wave file */
      PlayFile(NULL);
      break;
    case 7: /* Record wave file */
      open_window(frecord);
      break;
    case 8: /* Set mixer */
      printf("Not implemented\n");
      /* open_window(fmixer); */
      break;
    case 9: /* Messages */
      open_window(fmessages);
      break;
    case 10: /* Quit */;
      xv_end=1;
      break;
  }
}

typedef unsigned char DacPalette256[3][256];

void MakePalette(Xv_opaque obj)
{
#define PALETTE "anawave.xpl"
  FILE *pal;
  int i,j;
  DacPalette256 cores;

  if (COLOR_PALETTE) {
    pal=fopen(PALETTE,"rb");
    if (pal!=NULL) {
      fread(cores,sizeof(cores),1,pal);
      fclose(pal);
    }
    else printf("Palette file %s not found\n",PALETTE);
  }
  else
    for (i=0; i<256; i++)
      cores[0][i]=cores[1][i]=cores[2][i]=i;
  for (i=first_color;i<=last_color;i++) {
    j=255*(i-first_color)/(last_color-first_color);
    setrgbpalette(i,cores[0][j]>>2,cores[1][j]>>2,cores[2][j]>>2);
  }
}

/*
void MakePalette(Xv_opaque obj)
{
  s_int i,l0,l1,l2,l3,l4,l5,l6;
  DacPalette256 Hallo;

  l0=first_color;
  l6=last_color;
  if (COLOR_PALETTE) {
    l1=(7*l0+l6)/8;
    l2=(5*l0+3*l6)/8;
    l3=(l0+l6)/2;
    l4=(2*l0+6*l6)/8;
    l5=(l0+7*l6)/8;
    for(i=l0; i<=l6; i++) {
      if(i<=l1) {
        Hallo[i][0]=0;
        Hallo[i][1]=0;
        Hallo[i][2]=63*(i-l0)/(l1-l0);
      }
      if((i>l1)&&(i<=l2)) {
        Hallo[i][0]=0;
        Hallo[i][1]=63*(i-l1)/(l2-l1);
        Hallo[i][2]=63*(l2-i)/(l2-l1);
      }
      if((i>l2)&&(i<=l3)) {
        Hallo[i][0]=63*(i-l2)/(l3-l2);
        Hallo[i][1]=63;
        Hallo[i][2]=0;
      }
      if((i>l3)&&(i<=l4)) {
        Hallo[i][0]=63;
        Hallo[i][1]=63*(l4-i)/(l4-l3);
        Hallo[i][2]=63*(i-l3)/(l4-l3);
      }
      if((i>l4)&&(i<=l5)) {
        Hallo[i][0]=63*(l6-i)/(l6-l4);
        Hallo[i][1]=0;
        Hallo[i][2]=63*(l5-i)/(l5-l4);
      }
      if (i>l5) {
        Hallo[i][0]=63*(i-l4)/(l6-l4);
        Hallo[i][1]=0;
        Hallo[i][2]=0;
      }
    }
  }
  else
    for (i=l0; i<=l6; i++)
      Hallo[i][0]=Hallo[i][1]=Hallo[i][2]=63*(i-l0)/(l6-l0);
  for (i=first_color;i<=last_color;i++)
    setrgbpalette(i,Hallo[i][0],Hallo[i][1],Hallo[i][2]);
}
*/

s_int total,nomes;

void ReadDirectory(Xv_opaque obj)
{
  struct ffblk srec;
  short done;

  settextstyle(SMALL_FONT,HORIZ_DIR,4);
  if (obj!=cdirectory) {
    setfillstyle(SOLID_FILL,cdirectory->back_color);
    bar(0,0,cdirectory->dx,cdirectory->dy);
  }
  nomes=cdirectory->dx/78;
  total=0;
  done=findfirst(tmask->v.stextfield.panel_value,&srec,FA_DIREC|FA_ARCH|FA_RDONLY);
  while (!done) {
    outtextxy(total%nomes*78+3,total/nomes*8,srec.ff_name);
    done=findnext(&srec);
    total++;
  }
  total--;
}

void SelectFile(Xv_opaque obj)
{
  struct ffblk srec;
  char drive[MAXDRIVE];
  char dir[MAXDIR];
  char file[MAXFILE];
  char ext[MAXEXT];
  short i,k,done;

  if (ie_code!=MS_LEFT)
    return;
  k=(ie_locx-3)/78;
  if (k>=nomes)
    return;
  k+=(ie_locy-3)/8*nomes;
  if (k>total)
    return;
  i=0;
  done=findfirst(tmask->v.stextfield.panel_value,&srec,FA_DIREC|FA_ARCH|FA_RDONLY);
  while ((!done)&&(i<k)) {
    done=findnext(&srec);
    i++;
  }
  fnsplit(tmask->v.stextfield.panel_value,drive,dir,file,ext);
  sprintf(tfilename->v.stextfield.panel_value,"%s%s%s",drive,dir,srec.ff_name);
  update(tfilename);
}

void ApplyFilter(Xv_opaque obj)
{
  FILE *filt;
  double A[BMAX],C[BMAX],F[BMAX],L[BMAX],B[BMAX],H[BMAX],V1[BMAX],V2[BMAX],
         Vin,S,I2;
  s_int n,i,j,limit;
  char tst[3];
  unsigned l_int s,s1,s2;
  s_int sample;
  header h;
  double sr;

  if (!valid_file) {
    printf("Read the input file first, for verification\n");
    return;
  }
  if (!strncasecmp(get_value(tinfile),get_value(toutfile),12)) {
    printf("The input and output files must be different\n");
    return;
  }
  if ((filt=fopen(get_value(tfilter),"r"))==NULL) {
    printf("File %s not found\n",get_value(tfilter));
    return;
  }
  fscanf(filt,"%2s",tst);
  if (strcmp("SR",tst)) {
    printf("Invalid filter file\n");
    fclose(filt);
    return;
  }
  fscanf(filt,"%*c%lg",&sr);
  printf("Filter sampling rate: %lg Hz\n",sr);
  printf("Sample sampling rate: %lg Hz\n",sampling_rate);
  printf("Reading filter coefficients\n");
  fscanf(filt,"%*s%hd",&n);
  printf("The filter contains %d biquads\n",n);
  if (n>=BMAX) {
    printf("Too many biquads (max=%d)\n",BMAX);
    fclose(filt);
    return;
  }
  printf("Biquad coefficients:\n\n");
  for (i=0; i<n; i++) {
    fscanf(filt,"%*s%lg%*s%lg%*s%lg%*s%lg%*s%lg%*s%lg",&A[i],&C[i],&F[i],&L[i],&B[i],&H[i]);
    j=i+1;
    printf("A%d=%g\nC%d=%g\nF%d=%g\nL%d=%g\nB%d=%g\nH%d=%g\n\n",j,A[i],j,C[i],j,F[i],j,L[i],j,B[i],j,H[i]);
  }
  fclose(filt);
  if ((wave=fopen(get_value(tinfile),"rb"))==NULL) {
    printf("File %s not found\n",get_value(tinfile));
    return;
  }
  if ((filt=fopen(get_value(toutfile),"wb"))==NULL) {
    printf("Impossible to open file %s\n",get_value(toutfile));
    fclose(wave);
    return;
  }
  s1=(unsigned l_int)get_real(tfxmin);
  s2=s1+(unsigned l_int)get_real(tfxdelta);
  strcpy(h.riff,"RIFF");
  h.filesize=sizeof(header)-8+(s2-s1+1)*sample_size; /* The first 8 bytes do not count */
  strcpy(h.wave,"WAVE");
  strcpy(h.fmt,"fmt ");
  h.fmtsize=16;
  h.wFormatTag=1;
  h.nChannels=1;
  h.nSamplesPerSec=(l_int)sampling_rate;
  h.nAvgBytesPerSec=(l_int)sampling_rate*sample_size;
  h.nBlockAlign=sample_size;
  h.wBitsPerSample=sample_size*8;
  strcpy(h.data,"data");
  h.datasize=(s2-s1+1)*sample_size;
  fwrite(&h,sizeof(header),1,filt);
  fseek(wave,wave_start+(s1-1)*sample_size,0);
  printf("Filtering %s with filter %s...\n",get_value(tinfile),get_value(tfilter));
  if (sample_size==1) limit=127; else limit=32767;
  for (i=0; i<n; i++) V1[i]=V2[i]=0;
  for (s=s1; s<=s2; s++) {
    fread(&sample,sample_size,1,wave);
    if (sample_size==1) Vin=(sample & 255)-128;
    else Vin=sample;
    for (i=0; i<n; i++) {
      S=Vin+V2[i]+F[i]*V1[i];
      Vin=L[i]*V2[i]+B[i]*V1[i]-H[i]*S;
      I2=C[i]*V1[i];
      V1[i]-=A[i]*S;
      V2[i]+=I2;
    }
    if (Vin>limit) Vin=limit;
    else if (Vin<-limit) Vin=-limit;
    if (sample_size==1) sample=Vin+128; else sample=Vin;
    fwrite(&sample,sample_size,1,filt);
  }
  fclose(filt);
  fclose(wave);
  printf("Filtered sample saved as %s\n",get_value(toutfile));
}

void NotifyIn(Xv_opaque obj)
{
  set_value(tinfile,get_value(tfilename));
  update(tinfile);
  printf("Read the file to update this field\n");
}

void NotifyRecord(Xv_opaque obj)
{
  /* Notify handler for brecord */
  char *blaster, path[100], t1[100];
  s_int result;

  if ((blaster=getenv("SOUND"))==NULL) {
    printf("SOUND environment variable not found\n");
    return;
  }
  sprintf(path,"%s/%s",blaster,"record");
  sprintf(t1,"/A:%s /M:%s /R:%s /S:%s /Q",
    get_string(ssource),
    get_string(smode),
    get_string(ssize),
    get_string(ssampling));
  xv_set(mrecord,"Recording (Esc to quit)");
  result=spawnl(P_WAIT,path,"1",get_value(tsave),t1,NULL);
  if (result) printf("Some error occurred\n");
  xv_set(mrecord,"*");
  printf("Free memory: %lu\n",coreleft());
}

void CopyWLimits(Xv_opaque obj)
{
  set_real(tsdelta,get_real(txdelta));
  set_real(tsmin,get_real(txmin));
  update(tsmin);
  update(tsdelta);
}

void CopySLimits(Xv_opaque obj)
{
  set_real(txdelta,get_real(tsdelta));
  set_real(txmin,get_real(tsmin));
  update(txmin);
  update(txdelta);
}

void main(int argc,char *argv[])
{
  /* Inicialization */
  normal_bsize=10000;
  ComputeSinCos();
  xv_init(0,0);
  menu2=xv_create(NULL,MENU,
    XV_LABEL,"Options",
    ITEM_MENU,
      "Read signal",
      "Plot waveform",
      "Plot spectrogram",
      "Plot spectrum",
      "Filter signal",
      "Play signal",
      "Record signal",
      "Set mixer",
      "Messages",
      "Quit",
    NULL,
    NOTIFY_HANDLER,ProcessMenu2,
    NULL);
  /* Interface objects creation */
  fspectrum=xv_create(NULL,FRAME,
    XV_LABEL,"Spectrogram",
    Y,240,
    DX,639,
    DY,239,
    MENU_NAME,menu2,
    NULL);
  cspectrum=xv_create(fspectrum,CANVAS,
    BACK_COLOR,BLACK,
    FORE_COLOR,WHITE,
    MENU_NAME,menu2,
    NOTIFY_HANDLER,NotifySpectrum,
    EVENT_HANDLER,EventSpectrum,
    NULL);

  fpspectrum=xv_create(NULL,FRAME,
    XV_LABEL,"Spectrogram parameters",
    X,266,
    Y,23,
    DX,320,
    DY,204,
    MENU_NAME,menu2,
    NULL);
  tfftsize=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"FFT length",
    FIELD_TYPE,int_field,
    PANEL_INT,256,
    MIN_VALUE,1,
    NOTIFY_HANDLER,NFFTLength,
    NULL);
  tshift=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Shift step",
    Y,15,
    FIELD_TYPE,int_field,
    PANEL_INT,128,
    MIN_VALUE,1,
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  slog=xv_create(fpspectrum,SETTING,
    XV_LABEL,"Amplitude scale",
    Y,30,
    ITEM_SETTING,
      "Linear",
      "Sqrt",
      "Log",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,3,
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  swindow=xv_create(fpspectrum,SETTING,
    XV_LABEL,"Window",
    Y,45,
    ITEM_SETTING,
      "Rectang","Triang","Parab","Raised Cos",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,4,
    NOTIFY_HANDLER,ChangeWindow,
    NULL);
  tfmin=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Minimum frequency",
    Y,60,
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  tfmax=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Maximum frequency",
    Y,75,
    FIELD_TYPE,real_field,
    PANEL_REAL,"8000",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  sunit=xv_create(fpspectrum,SETTING,
    XV_LABEL,"Frequency unit",
    Y,90,
    ITEM_SETTING,
      "Hertz",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  spalette=xv_create(fpspectrum,SETTING,
    XV_LABEL,"Palette",
    Y,90,
    X,170,
    ITEM_SETTING,
      "Color",
      "B&W",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NOTIFY_HANDLER,MakePalette,
    NULL);
  tsmin=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"First sample ",
    Y,105,
    FIELD_TYPE,real_field,
    PANEL_REAL,"1",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  tsdelta=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Segments",
    Y,120,
    FIELD_TYPE,real_field,
    PANEL_REAL,"9999",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  tamin=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Minimum amplitude",
    Y,135,
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  tamax=xv_create(fpspectrum,TEXTFIELD,
    XV_LABEL,"Maximum amplitude",
    Y,150,
    FIELD_TYPE,real_field,
    PANEL_REAL,"40000",
    NOTIFY_HANDLER,InvalidateSpectrum,
    NULL);
  bapplyspectrum=xv_create(fpspectrum,BUTTON,
    XV_LABEL,"Apply",
    Y,165,
    NOTIFY_HANDLER,PlotSpectrum,
    NULL);
  xv_create(fpspectrum,BUTTON,
    XV_LABEL,"Copy waveform limits",
    Y,165, X,50,
    NOTIFY_HANDLER,CopyWLimits,
    NULL);

  fwave=xv_create(NULL,FRAME,
    XV_LABEL,"Waveform",
    DX,639,
    DY,239,
    MENU_NAME,menu2,
    NULL);
  cwave=xv_create(fwave,CANVAS,
    BACK_COLOR,BLACK,
    FORE_COLOR,WHITE,
    NOTIFY_HANDLER,NotifyWave,
    EVENT_HANDLER,EventWave,
    MENU_NAME,menu2,
    NULL);
  fwave->v.sframe.mouse_obj=fwave;
  fpwave=xv_create(NULL,FRAME,
    XV_LABEL,"Waveform Parameters",
    X,6,
    Y,23,
    DX,259,
    DY,115,
    MENU_NAME,menu2,
    NULL);
  tymax=xv_create(fpwave,TEXTFIELD,
    XV_LABEL,"Maximum signal",
    FIELD_TYPE,real_field,
    PANEL_REAL,"32767",
    NOTIFY_HANDLER,InvalidateWave,
    NULL);
  tymin=xv_create(fpwave,TEXTFIELD,
    XV_LABEL,"Minimum signal",
    Y,15,
    FIELD_TYPE,real_field,
    PANEL_REAL,"-32767",
    NOTIFY_HANDLER,InvalidateWave,
    NULL);
  txmin=xv_create(fpwave,TEXTFIELD,
    XV_LABEL,"First sample ",
    Y,30,
    FIELD_TYPE,real_field,
    NOTIFY_HANDLER,InvalidateWave,
    NULL);
  txdelta=xv_create(fpwave,TEXTFIELD,
    XV_LABEL,"Segments",
    Y,45,
    FIELD_TYPE,real_field,
    PANEL_REAL,"10000",
    NOTIFY_HANDLER,InvalidateWave,
    NULL);
  stime=xv_create(fpwave,SETTING,
    XV_LABEL,"Time unit",
    Y,60,
    ITEM_SETTING,"samples","seconds",NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NOTIFY_HANDLER,InvalidateWave,
    NULL);
  bapplywave=xv_create(fpwave,BUTTON,
    XV_LABEL,"Apply",
    Y,75,
    NOTIFY_HANDLER,PlotWave,
    NULL);
  xv_create(fpwave,BUTTON,
    XV_LABEL,"Copy spectrogram limits",
    Y,75, X,50,
    NOTIFY_HANDLER,CopySLimits,
    NULL);

  fio=xv_create(NULL,FRAME,
    XV_LABEL,"Files",
    DX,319,
    DY,120,
    MENU_NAME,menu2,
    NULL);
  tfilename=xv_create(fio,TEXTFIELD,
    XV_LABEL,"File name",
    VALUE_LENGTH,28,
    PANEL_VALUE,"x.wav",
    NULL);
  bread=xv_create(fio,BUTTON,
    XV_LABEL,"Read",
    Y,15,
    NOTIFY_HANDLER,ReadFile,
    NULL);
  bwrite=xv_create(fio,BUTTON,
    XV_LABEL,"Record",
    X,41,
    Y,15,
    NOTIFY_HANDLER,WriteFile,
    NULL);
   bplay=xv_create(fio,BUTTON,
    XV_LABEL,"Play",
    X,98,
    Y,15,
    NOTIFY_HANDLER,PlayFile,
    NULL);
  tmask=xv_create(fio,TEXTFIELD,
    XV_LABEL,"Mask",
    PANEL_VALUE,"*.wav",
    VALUE_LENGTH,33,
    NOTIFY_HANDLER,ReadDirectory,
    Y,30,
    NULL);
  cdirectory=xv_create(fio,CANVAS,
    Y,45,
    NOTIFY_HANDLER,ReadDirectory,
    EVENT_HANDLER,SelectFile,
    NULL);
  fio->v.sframe.mouse_obj=bread;

  fmessages=xv_create(NULL,FRAME,
    XV_LABEL,"Messages",
    DX,319,
    DY,239,
    X,320,
    MENU_NAME,menu2,
    NULL);
  tty1=xv_create(fmessages,TTY,
    NULL);

  fpfft=xv_create(NULL,FRAME,
    XV_LABEL,"Spectrum parameters",
    X,266,
    Y,23,
    DX,290,
    DY,114,
    MENU_NAME,menu2,
    NULL);
  tfftfmin=xv_create(fpfft,TEXTFIELD,
    XV_LABEL,"Minimum frequency",
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NOTIFY_HANDLER,InvalidateFFT,
    NULL);
  tfftfmax=xv_create(fpfft,TEXTFIELD,
    XV_LABEL,"Maximum frequency",
    Y,15,
    FIELD_TYPE,real_field,
    PANEL_REAL,"8000",
    NOTIFY_HANDLER,InvalidateFFT,
    NULL);
  tfftsmin=xv_create(fpfft,TEXTFIELD,
    XV_LABEL,"First sample ",
    Y,30,
    FIELD_TYPE,real_field,
    PANEL_REAL,"1",
    NOTIFY_HANDLER,InvalidateFFT,
    NULL);
  tfftamin=xv_create(fpfft,TEXTFIELD,
    XV_LABEL,"Minimum amplitude",
    Y,45,
    FIELD_TYPE,real_field,
    PANEL_REAL,"0",
    NOTIFY_HANDLER,InvalidateFFT,
    NULL);
  tfftamax=xv_create(fpfft,TEXTFIELD,
    XV_LABEL,"Maximum amplitude",
    Y,60,
    FIELD_TYPE,real_field,
    PANEL_REAL,"40000",
    NOTIFY_HANDLER,InvalidateFFT,
    NULL);
  bapplyfft=xv_create(fpfft,BUTTON,
    XV_LABEL,"Plot",
    Y,75,
    NOTIFY_HANDLER,PlotFFT,
    NULL);
  blistfft=xv_create(fpfft,BUTTON,
    XV_LABEL,"List",
    Y,75,X,40,
    NOTIFY_HANDLER,ListFFT,
    NULL);

  ffft=xv_create(NULL,FRAME,
    XV_LABEL,"Spectrum",
    DX,319,
    DY,239,
    MENU_NAME,menu2,
    NULL);
  cfft=xv_create(ffft,CANVAS,
    NOTIFY_HANDLER,NotifyFFT,
    EVENT_HANDLER,EventFFT,
    MENU_NAME,menu2,
    NULL);

   ffilter=xv_create(NULL,FRAME,
    XV_LABEL,"Filter",
    DX,285,
    DY,114,
    MENU_NAME,menu2,
    NULL);
  tinfile=xv_create(ffilter,TEXTFIELD,
    XV_LABEL,"Input wave file  ",
    PANEL_VALUE,"filter0.mul",
    NOTIFY_HANDLER,NotifyIn,
    NULL);
  tfilter=xv_create(ffilter,TEXTFIELD,
    XV_LABEL,"Filter mult. file",
    PANEL_VALUE,"filter1.mul",
    Y,15,
    NULL);
  toutfile=xv_create(ffilter,TEXTFIELD,
    XV_LABEL,"Output wave file ",
    PANEL_VALUE,"filter1.wav",
    Y,30,
    NULL);
  tfxmin=xv_create(ffilter,TEXTFIELD,
    XV_LABEL,"First sample ",
    Y,45,
    FIELD_TYPE,real_field,
    NULL);
  tfxdelta=xv_create(ffilter,TEXTFIELD,
    XV_LABEL,"Segments",
    Y,60,
    FIELD_TYPE,real_field,
    NULL);
  bapplyfilter=xv_create(ffilter,BUTTON,
    XV_LABEL,"Apply",
    Y,75,
    NOTIFY_HANDLER,ApplyFilter,
    NULL);

  frecord=xv_create(NULL,FRAME,
    XV_LABEL,"Sound Recording",
    X,41,
    Y,25,
    DX,262,
    DY,115,
    NULL);
  ssource=xv_create(frecord,SETTING,
    XV_LABEL,"Source",
    Y,15,
    ITEM_SETTING,
      "MIC",
      "CD",
      "LINE",
      "FM",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NULL);
  smode=xv_create(frecord,SETTING,
    XV_LABEL,"Mode",
    Y,30,
    ITEM_SETTING,
      "MONO",
      "STEREO",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NULL);
  ssampling=xv_create(frecord,SETTING,
    XV_LABEL,"Sampling rate",
    Y,60,
    ITEM_SETTING,
      "11025",
      "22050",
      "44100",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,1,
    NULL);
  tsave=xv_create(frecord,TEXTFIELD,
    XV_LABEL,"Save as (.wav)",
    PANEL_VALUE,"out.wav",
    NULL);
  brecord=xv_create(frecord,BUTTON,
    XV_LABEL,"Start",
    Y,75,
    NOTIFY_HANDLER,NotifyRecord,
    NULL);
  mrecord=xv_create(frecord,MESSAGE,
    XV_LABEL,"*",
    FORE_COLOR,4,
    X,52,
    Y,75,
    NULL);
  ssize=xv_create(frecord,SETTING,
    XV_LABEL,"Sample size (bits)",
    Y,45,
    ITEM_SETTING,
      "8",
      "16",
    NULL,
    EXCLUSIVE,1,
    SEL_SETTING,2,
    NULL);

  fmixer=xv_create(NULL,FRAME,
    XV_LABEL,"Mixer levels",
    DX,245,
    DY,69,
    DYMIN,27,
    NULL);
  tmaster=xv_create(fmixer,TEXTFIELD,
    XV_LABEL,"Master    ",
    VALUE_LENGTH,5,
    FIELD_TYPE,int_field,
    PANEL_INT,127,
    MIN_VALUE,0,
    MAX_VALUE,255,
    NULL);
  tmic=xv_create(fmixer,TEXTFIELD,
    XV_LABEL,"Microphone",
    Y,14,
    VALUE_LENGTH,5,
    FIELD_TYPE,int_field,
    PANEL_INT,200,
    MIN_VALUE,0,
    MAX_VALUE,255,
    NULL);
  tbalance=xv_create(fmixer,TEXTFIELD,
    XV_LABEL,"Balance",
    X,137,
    VALUE_LENGTH,4,
    FIELD_TYPE,int_field,
    MIN_VALUE,-5,
    MAX_VALUE,5,
    NULL);
  sagc=xv_create(fmixer,SETTING,
    XV_LABEL,"AGC",
    X,137,
    Y,15,
    ITEM_SETTING,
      " ",
    NULL,
    NULL);
  bapplymixer=xv_create(fmixer,BUTTON,
    XV_LABEL,"Apply",
    Y,30,
    NOTIFY_HANDLER,NotifyMixer,
    NULL);

  open_window(fspectrum);
  open_window(fwave);
  printf("AnaWave, Version %s\nBy Antonio Carlos M. de Queiroz\nCOPPE/UFRJ\nE-mail: acmq@coe.ufrj.br\n",version);
  MakePalette(NULL);
  xv_main_loop(fio);
  /* Exit */
  /* restorecrtmode(); */
}
