C
C     An example FORTRAN program using the eps09s.cpp
C
      PROGRAM main
        IMPLICIT NONE
        INTEGER i, order, set, A
        DOUBLE PRECISION x, q, s, ta, ta_eps09s, r_g, r_uv, r_us
        DOUBLE PRECISION ruv, rdv, ru, rd, rs, rc, rb, rg
        DOUBLE PRECISION cuv(4), cdv(4), cu(4), cd(4), cs(4), cc(4), 
     *    cb(4), cg(4)
        order = 1
        set = 1
        A = 208
        x = 0.02
        q = 1.4
        s = 5
        CALL eps09s(order,set,A,x,Q,s,ruv,rdv,ru,rd,rs,rc,rb,rg)
        ta = ta_eps09s(A,s)
        PRINT*, rg/ta, ruv/ta, ru/ta
        r_g = ta
        r_uv = ta
        r_us = ta
        CALL eps09s_c(order,set,A,x,Q,cuv,cdv,cu,cd,cs,cc,cb,cg)
        DO i = 1,4
           r_g = r_g + cg(i)*ta**(i+1)
           r_uv = r_uv + cuv(i)*ta**(i+1)
           r_us = r_us + cu(i)*ta**(i+1)
        END DO
        PRINT*, r_g/ta, r_uv/ta, r_us/ta

        CALL test_eps09s()

      END PROGRAM
