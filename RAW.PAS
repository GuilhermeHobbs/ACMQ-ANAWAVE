{
Este programa escreve um arquivo .raw de som.
Nao ha cabecalho.
Em 8 bits devem ser escritos bytes, com 128 sendo o centro (0).
Em 16 bits, basta escrever inteiros.
Em estereo, parece que os canais se alternam.
}

CONST
  inicio=10000;
  final=40000;
  f1=5000;
  f2=0;
  a1=20000;
  amostragem=11025;
VAR
  saida:FILE of INTEGER;
  i:LONGINT;
  v:INTEGER;
  dpi,t,wa:EXTENDED;
  a,b:EXTENDED;
BEGIN
  dpi:=2.0*Pi;
  Assign(saida,'teste.raw');
  ReWrite(saida);
  a:=dpi*(f2-f1)/(final-inicio);
  b:=dpi*f1-a*inicio;
  FOR i:=inicio TO final DO BEGIN
    wa:=a/2*i+b; {/2!}
    t:=i/amostragem;
    v:=Round(a1*cos(wa*t));
    Write(saida,v);
  END;
  Close(saida);
END.