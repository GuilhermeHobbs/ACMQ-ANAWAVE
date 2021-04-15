VAR
  A,B:ARRAY[0..255] OF REAL;
  i,n,m:INTEGER;
  DoisPi:REAL;

FUNCTION Inverso(x:INTEGER):INTEGER;
VAR
  i,u:INTEGER;
BEGIN
  u:=0;
  i:=n shr 1;
  REPEAT
    IF Odd(x) THEN u:=u+i;
    i:=i shr 1;
    x:=x shr 1
  UNTIL x=0;
  Inverso:=u
END; {Inverso}

PROCEDURE FFT(n,m:INTEGER; direto:BOOLEAN);  {"Fast Fourier Transform"}
VAR
  k,k1,l,s,i,j,i1,i2:INTEGER;
  x1,x2,y1,y2:REAL;
BEGIN
  FOR k:=m-1 DOWNTO 0 DO BEGIN
    k1:=1 shl k;
    l:=0;
    REPEAT
      s:=Inverso(l div k1);
      x1:=Cos(DoisPi*s/n); y1:=Sin(DoisPi*s/n);
      IF direto THEN y1:=-y1;
      FOR j:=0 TO k1-1 DO BEGIN
        i1:=j+l+k1; i2:=j+l;
        x2:=A[i1]*x1-B[i1]*y1;
        y2:=A[i1]*y1+B[i1]*x1;
        A[i2]:=a[i2]+x2;
        B[i2]:=B[i2]+y2;
        A[i1]:=A[i2]-x2-x2;
        B[i1]:=B[i2]-y2-y2;
      END;
      l:=l+k1 shl 1
    UNTIL l>n-1
  END;
  FOR i:=0 TO n-1 DO BEGIN
    s:=Inverso(i);
    IF s>i THEN BEGIN
      x1:=A[i]; A[i]:=A[s]; A[s]:=x1;
      x1:=B[i]; B[i]:=B[s]; B[s]:=x1;
    END
  END;
  IF direto THEN
    FOR i:=0 TO n-1 DO BEGIN A[i]:=A[i]/n; B[i]:=B[i]/n END;
END; {FFT}

BEGIN
  DoisPi:=2*Pi;
  Write('Ordem da FFT: '); ReadLn(m);
  n:=1 shl m;
  WriteLn('FFT com ',n,' pontos');
  FOR i:=0 TO n-1 DO BEGIN
    Write('A(',i,') B(',i,'): '); ReadLn(A[i],B[i]);
  END;
  WriteLn;
  WriteLn('FFT direta:');
  FFT(n,m,TRUE);
  FOR i:=0 TO n-1 DO WriteLn('A(',i,') B(',i,'): ',A[i]:10:5,B[i]:10:5,'j');
  ReadLn;
  WriteLn('FFT inversa:');
  FFT(n,m,FALSE);
  FOR i:=0 TO n-1 DO WriteLn('A(',i,') B(',i,'): ',A[i]:10:5,B[i]:10:5,'j');
  ReadLn;
END.