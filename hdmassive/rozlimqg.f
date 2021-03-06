      FUNCTION HQHDqgM0(T1,RO,NL)
      IMPLICIT DOUBLE PRECISION (A-Z)
      INTEGER NL
      DATA PI/3.141 592 653 589 793/
      T2 = 1-T1
      XLRO = LOG(RO/4.0D0)
      DD = 0
      DD = DD
      HQHDqgM0 = DD
      RETURN
      END
      FUNCTION ASHPqgM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      XLRO = LOG(RO/4.0D0)
      FOUR = 4
      XLOG4 = LOG(FOUR)
      T2 = -TX-T1+1
      PP = 2*XLRO+2*XLOG4+LOG(((-T2-T1+1)**2+2*T1*(-T2-T1+1)+T1**2)/((-T
     1   2-T1+1)**2+(2*T1-2)*(-T2-T1+1)+T1**2-2*T1+1)/1.6D1)
      PP = (-5.0D0)*PP*(T2+T1-1)*(3*T2**2+6*T1*T2-4*T2+4*T1**2-4*T1+2)/(
     1   7.2D1*(T2-1)*T2)
      TMP0 = XLRO+XLOG4-LOG(-4*(T1-1)*(1-T2))
      TMP0 = (-5.0D0)*(-T2-T1+1)*(2*T2**2+4*T1*T2-4*T2+4*T1**2-6*T1+3)*T
     1   MP0/(3.6D1*(T2-1)*T2)
      PP = TMP0+PP
      TMP0 = XLRO+XLOG4-LOG(4*((-T2-T1+1)**2-2*(-T2-T1+1)+1))
      TMP0 = 5.0D0*(-T2-T1+1)*(T2**2+3*T1*T2-2*T2+2*T1**2-T1)*TMP0/(7.2D
     1   1*(-T2-T1)*T2)
      PP = TMP0+PP
      PP = 5.0D0*(-T2-T1+1)*(T2+T1)**2*(XLRO+XLOG4-LOG(4*((-T2-T1+1)**2-
     1   2*(-T2-T1+1)+1)))/(7.2D1*(-T2-T1)**3)+PP
      PP = (-5.0D0)*(-T2-T1+1)**2*XLRO/(7.2D1*(T1-1)*T2)+PP
      PP = 5.0D0*(-T2-T1+1)*(T2**2+3*T1*T2+3*T1-2)*XLRO/(3.6D1*(1-T2)*T2
     1   )+PP
      PP = 5.0D0*(-T2-T1+1)*(T2-1)*(T1*T2+T1-1)*XLRO/(3.6D1*(1-T2)**2*T2
     1   )+PP
      PP = 5.0D0*(-T2-T1+1)*(T1*T2-2*T2+2*T1**2-4*T1+2)*XLRO/(7.2D1*(T1-
     1   1)*T2)+PP
      PP = 5.0D0*(-T2-T1+1)*(T2**2+2*T1*T2-2*T2+4*T1**2-4*T1+2)*LOG(-T2/
     1   (T1-1))/(3.6D1*(T2-1)*T2)+PP
      PP = (-5.0D0)*(-T2-T1+1)*(T2+T1)*(T1*T2+T2+T1**2-T1)/(7.2D1*(-T2-T
     1   1)**2*T2)+PP
      PP = (-5.0D0)*(-T2-T1+1)*(T2+2*T1-1)/(3.6D1*(1-T2)*T2)+PP
      PP = 5.0D0*LOG(T1/(1-T2))*(-T2-T1+1)*(3*T2**2+6*T1*T2-4*T2+8*T1**2
     1   -6*T1+3)/(3.6D1*(T2-1)*T2)+PP
      PP = 5.0D0*(T2+T1-1)*(T1*T2+4*T1**2-6*T1+2)/(1.44D2*(T1-1)*T2)+PP
      PP = (-5.0D0)*(T2+T1-1)*(T1*T2**3-4*T2**3+6*T1**2*T2**2-8*T1*T2**2
     1   +8*T2**2-4*T1**2*T2+13*T1*T2-12*T2+14*T1**2-22*T1+8)/(1.44D2*(T
     2   1-1)*(T2-1)**2*T2)+PP
      ASHPqgM0 = PP
      RETURN 
      END 
      FUNCTION HQHLqgM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      T2 = -TX-T1+1
      T11 = 1/(1-T1)
      T22 = 1/(TX+T1)
      PP = (9*T1+7)*(-T2-T1+1)**2/((T1-1)**2*T2)/7.2D1
      PP = (-T2-T1+1)*(9*T2**2+23*T1*T2-16*T2-18*T1**2+27*T1)/((-T2-T1)*
     1   T2)/7.2D1+PP
      PP = (-T2-T1+1)*(T2+T1)*(7*T2-9*T1)/(-T2-T1)**3/7.2D1+PP
      PP = (-T2-T1+1)*(7*T2**2+23*T1*T2+2*T2+23*T1-16)/((1-T2)*T2)/3.6D1
     1   +PP
      PP = (-T2-T1+1)*(T2-1)*(9*T1*T2+2*T2+9*T1-9)/((1-T2)**2*T2)/3.6D1+
     1   PP
      PP = PP-(-T2-T1+1)*(27*T1*T2+16*T2+18*T1**2+14*T1-16)/((T1-1)*T2)/
     1   7.2D1
      PP = PP-(T2-2)*(T2+T1-1)*(3*T2**2+6*T1*T2-4*T2+4*T1**2-4*T1+2)/((T
     1   2-1)*T2**2)/2.0D0
      TMP0 = 16*T1**2*T2**4*T22**2-16*T1**2*T2**3*T22**2+8*T1**2*T2**2*T
     1   22**2-16*T1*T2**4*T22+16*T1*T2**3*T22-8*T1*T2**2*T22+8*T2**4-8*
     2   T2**3-18*T1**5*T11**3*T2**2+54*T1**4*T11**3*T2**2-71*T1**3*T11*
     3   *3*T2**2+52*T1**2*T11**3*T2**2-21*T1*T11**3*T2**2+4*T11**3*T2**
     4   2+4*T2**2+36*T1**5*T11**2*T2-108*T1**4*T11**2*T2+142*T1**3*T11*
     5   *2*T2-104*T1**2*T11**2*T2+42*T1*T11**2*T2-8*T11**2*T2+36*T1**6*
     6   T11**2-144*T1**5*T11**2+250*T1**4*T11**2-246*T1**3*T11**2+146*T
     7   1**2*T11**2-50*T1*T11**2+8*T11**2
      TMP0 = -(T2+T1-1)*TMP0/(T1*T2**2)/9.0D0
      PP = TMP0+PP
      HQHLqgM0 = PP
      RETURN 
      END 
      FUNCTION HQHPqgM0(TX,T1,RO)
      IMPLICIT DOUBLE PRECISION (A-Z)
      XLRO = LOG(RO/4.0D0)
      FOUR = 4
      XLOG4 = LOG(FOUR)
      T2 = -TX-T1+1
      T11 = 1/(1-T1)
      T22 = 1/(TX+T1)
      PP = 2*XLRO+2*XLOG4+LOG(((-T2-T1+1)**2+2*T1*(-T2-T1+1)+T1**2)/((-T
     1   2-T1+1)**2+(2*T1-2)*(-T2-T1+1)+T1**2-2*T1+1)/1.6D1)
      PP = PP*(T2-2)*(T2+T1-1)*(3*T2**2+6*T1*T2-4*T2+4*T1**2-4*T1+2)/((T
     1   2-1)*T2**2)/8.0D0
      TMP0 = XLRO+XLOG4-LOG(-4*(T1-1)*(1-T2))
      TMP0 = -(-T2-T1+1)*(2*T2-1)*(2*T2**2+4*T1*T2-4*T2+4*T1**2-6*T1+3)*
     1   TMP0/((T2-1)*T2)/4.0D0
      PP = TMP0+PP
      TMP0 = XLRO+XLOG4-LOG(4*((-T2-T1+1)**2-2*(-T2-T1+1)+1))
      TMP0 = -(-T2-T1+1)*(9*T2**2+23*T1*T2-16*T2-18*T1**2+27*T1)*TMP0/((
     1   -T2-T1)*T2)/7.2D1
      PP = TMP0+PP
      PP = PP-(-T2-T1+1)*(T2+T1)*(7*T2-9*T1)*(XLRO+XLOG4-LOG(4*((-T2-T1+
     1   1)**2-2*(-T2-T1+1)+1)))/(-T2-T1)**3/7.2D1
      PP = PP-(9*T1+7)*(-T2-T1+1)**2*XLRO/((T1-1)**2*T2)/7.2D1
      PP = PP-(-T2-T1+1)*(7*T2**2+23*T1*T2+2*T2+23*T1-16)*XLRO/((1-T2)*T
     1   2)/3.6D1
      PP = PP-(-T2-T1+1)*(T2-1)*(9*T1*T2+2*T2+9*T1-9)*XLRO/((1-T2)**2*T2
     1   )/3.6D1
      PP = (-T2-T1+1)*(27*T1*T2+16*T2+18*T1**2+14*T1-16)*XLRO/((T1-1)*T2
     1   )/7.2D1+PP
      TMP0 = 16*T1**2*T2**4*T22**2-16*T1**2*T2**3*T22**2+8*T1**2*T2**2*T
     1   22**2-16*T1*T2**4*T22+16*T1*T2**3*T22-8*T1*T2**2*T22+8*T2**4-8*
     2   T2**3-18*T1**5*T11**3*T2**2+54*T1**4*T11**3*T2**2-71*T1**3*T11*
     3   *3*T2**2+52*T1**2*T11**3*T2**2-21*T1*T11**3*T2**2+4*T11**3*T2**
     4   2+4*T2**2+36*T1**5*T11**2*T2-108*T1**4*T11**2*T2+142*T1**3*T11*
     5   *2*T2-104*T1**2*T11**2*T2+42*T1*T11**2*T2-8*T11**2*T2+36*T1**6*
     6   T11**2-144*T1**5*T11**2+250*T1**4*T11**2-246*T1**3*T11**2+146*T
     7   1**2*T11**2-50*T1*T11**2+8*T11**2
      TMP0 = (T2+T1-1)*TMP0*XLRO/(T1*T2**2)/9.0D0
      PP = TMP0+PP
      TMP0 = 36*T1*T2**4+144*T1**2*T2**3-117*T1*T2**3+16*T2**3+216*T1**3
     1   *T2**2-378*T1**2*T2**2+212*T1*T2**2-48*T2**2+144*T1**4*T2-468*T
     2   1**3*T2+496*T1**2*T2-258*T1*T2+64*T2-144*T1**4+288*T1**3-280*T1
     3   **2+136*T1-32
      TMP0 = (-T2-T1+1)*LOG(-T2/(T1-1))*TMP0/(T1*(T2-1)*T2**2)/3.6D1
      PP = TMP0+PP
      PP = (-T2-T1+1)*(32*T2**4+64*T1*T2**3-32*T2**3+32*T1**2*T2**2-25*T
     1   1*T2**2+7*T2**2-2*T1**2*T2-16*T1*T2-9*T1**3+9*T1**2)/((-T2-T1)*
     2   *2*T2)/7.2D1+PP
      PP = PP-(-T2-T1+1)*(2*T2**2+4*T1*T2+5*T2+14*T1-23)/((1-T2)*T2)/3.6
     1   D1
      TMP0 = 32*T2**5+100*T1*T2**4-64*T2**4+136*T1**2*T2**3-159*T1*T2**3
     1   +64*T2**3+72*T1**3*T2**2-154*T1**2*T2**2+136*T1*T2**2-64*T2**2-
     2   72*T1**3*T2+154*T1**2*T2-159*T1*T2+64*T2+72*T1**3-136*T1**2+100
     3   *T1-32
      TMP0 = LOG(T1/(1-T2))*(-T2-T1+1)*TMP0/(T1*(T2-1)*T2**2)/3.6D1
      PP = TMP0+PP
      PP = PP-(T2+T1-1)*(9*T1*T2+22*T2-36*T1**2+8*T1-4)/((T1-1)*T2)/1.44
     1   D2
      TMP0 = 64*T1**2*T2**6-128*T1*T2**6+64*T2**6+464*T1**3*T2**5-1120*T
     1   1**2*T2**5+848*T1*T2**5-192*T2**5+544*T1**4*T2**4-2239*T1**3*T2
     2   **4+3001*T1**2*T2**4-1434*T1*T2**4+160*T2**4+576*T1**5*T2**3-23
     3   74*T1**4*T2**3+4132*T1**3*T2**3-3436*T1**2*T2**3+1038*T1*T2**3+
     4   288*T1**6*T2**2-2016*T1**5*T2**2+4420*T1**4*T2**2-4579*T1**3*T2
     5   **2+2485*T1**2*T2**2-534*T1*T2**2-32*T2**2-576*T1**6*T2+2304*T1
     6   **5*T2-3342*T1**4*T2+2414*T1**3*T2-1138*T1**2*T2+338*T1*T2+288*
     7   T1**6-864*T1**5+992*T1**4-672*T1**3+384*T1**2-128*T1
      TMP0 = (T2+T1-1)*TMP0/((T1-1)**2*T1*(T2-1)**2*T2**2)/1.44D2
      PP = TMP0+PP
      HQHPqgM0 = PP
      RETURN 
      END 
