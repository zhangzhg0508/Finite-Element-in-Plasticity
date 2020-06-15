#include <fem.hpp> // Fortran EMulation library of fable module

namespace epca {

using namespace fem::major_types;

using fem::common;

//C$DEBUG
//C$LARGE
//C(Page 210)
void
algor(
  arr_ref<float> fixed,
  int const& iincs,
  int const& iiter,
  int& kresl,
  int const& mtotv,
  int const& nalgo,
  int const& ntotv)
{
  fixed(dimension(mtotv));
  //C******************************************************************
  //C
  //C*** THIS SUBROUTINE SETS EQUATION RESOLUTION INDEX,KREST
  //C
  //C******************************************************************
  kresl = 2;
  if (nalgo == 1 && iincs == 1 && iiter == 1) {
    kresl = 1;
  }
  if (nalgo == 2) {
    kresl = 1;
  }
  if (nalgo == 3 && iiter == 1) {
    kresl = 1;
  }
  if (nalgo == 4 && iincs == 1 && iiter == 1) {
    kresl = 1;
  }
  if (nalgo == 4 && iiter == 2) {
    kresl = 1;
  }
  if (iiter == 1) {
    return;
  }
  int itotv = fem::int0;
  FEM_DO_SAFE(itotv, 1, ntotv) {
    fixed(itotv) = 0.0f;
  }
}

//C
//C$DEBUG
//C$LARGE
//C(Page 201)
void
echo(
  common& cmn)
{
  common_read read(cmn);
  common_write write(cmn);
  arr_1d<1, int> ntitl(fem::fill0);
  static const char* format_900 =
    "(/,/,' NOW FOLLOWS A LISTING OF POST-DISASTER DATA CARDS',/)";
  static const char* format_910 = "(20x,i1)";
  //C*****************************************************************
  //C
  //C*** IF DATA ERRORS HAVE BEEN DETECTED BY SUBROUTINES CHECK1 OR
  //C    CHECK2,THIS SUBROUTINE READS AND WRITE THE REMAINING DATA CARDS
  //C
  //C********************************************************************
  write(6, format_900);
  //write(6, format_900);
  statement_10:
  read(5, "(i1)"), ntitl;
  write(6, format_910), ntitl;
  //write(6, format_910), ntitl;
  goto statement_10;
}

//C$DEBUG
//C$LARGE
//C(Page 200)
void
check1(
  common& cmn,
  int const& ndofn,
  int const& nelem,
  int const& ngaus,
  int const& nmats,
  int const& nnode,
  int const& npoin,
  int const& nstre,
  int const& ntype,
  int const& nvfix,
  int const& ncrit,
  int const& nalgo,
  int const& nincs)
{
  common_write write(cmn);
  int ieror = fem::int0;
  arr_1d<24, int> neror(fem::fill0);
  int keror = fem::int0;
  static const char* format_900 = "(/,/,' *** DIAGNOSIS BY CHECK1, ERROR',i3)";
  //C**********************************************************************
  //C
  //C*** THIS SUBROUTINE CHECKS THE MAIN CONTROL DATA
  //C
  //C**********************************************************************
  FEM_DO_SAFE(ieror, 1, 12) {
    neror(ieror) = 0;
  }
  //C
  //C*** CREATE THE DIAGNOSTIC MESSAGES
  //C
  if (npoin <= 0) {
    neror(1) = 1;
  }
  if (nelem * nnode < npoin) {
    neror(2) = 1;
  }
  if (nvfix < 2 || nvfix > npoin) {
    neror(3) = 1;
  }
  if (nincs < 1) {
    neror(4) = 1;
  }
  if (ntype < 1 || ntype > 3) {
    neror(5) = 1;
  }
  if (nnode < 4 || nnode > 9) {
    neror(6) = 1;
  }
  if (ndofn < 2 || ndofn > 5) {
    neror(7) = 1;
  }
  if (nmats < 1 || nmats > nelem) {
    neror(8) = 1;
  }
  if (ncrit < 1 || ncrit > 4) {
    neror(9) = 1;
  }
  if (ngaus < 2 || ngaus > 3) {
    neror(10) = 1;
  }
  if (nalgo < 1 || nalgo > 4) {
    neror(11) = 1;
  }
  if (nstre < 3 || nstre > 5) {
    neror(12) = 1;
  }
  //C
  //C*** EITHER RETURN, OR ELSE PRINT THE ERRORS DIAGNOSED
  //C
  keror = 0;
  FEM_DO_SAFE(ieror, 1, 12) {
    if (neror(ieror) == 0) {
      goto statement_20;
    }
    keror = 1;
    write(6, format_900), ieror;
    //write(6, format_900), ieror;
    statement_20:;
  }
  if (keror == 0) {
    return;
  }
  //C
  //C*** OTHERWISE ECHO ALL THE REMAINIGN DATA WITHOUT FURTHER COMMEN
  //C
  echo(cmn);
}

//C
//C$DEBUG
//C$LARGE
//C(Page 202)
void
check2(
  common& cmn,
  arr_cref<float, 2> coord,
  arr_cref<int> iffix,
  arr_ref<int, 2> lnods,
  arr_cref<int> matno,
  int const& melem,
  int const& mfron,
  int const& mpoin,
  int const& mtotv,
  int const& mvfix,
  arr_ref<int> ndfro,
  int const& ndofn,
  int const& nelem,
  int const& nmats,
  int const& nnode,
  arr_cref<int> nofix,
  int const& npoin,
  int const& nvfix)
{
  coord(dimension(mpoin, 2));
  iffix(dimension(mtotv));
  lnods(dimension(melem, 9));
  matno(dimension(melem));
  ndfro(dimension(melem));
  nofix(dimension(mvfix));
  common_write write(cmn);
  int ieror = fem::int0;
  arr_1d<24, int> neror(fem::fill0);
  int ielem = fem::int0;
  int ipoin = fem::int0;
  int kpoin = fem::int0;
  int jpoin = fem::int0;
  int idime = fem::int0;
  int inode = fem::int0;
  int kstar = fem::int0;
  int kzero = fem::int0;
  int klast = fem::int0;
  int nlast = fem::int0;
  float sigma = fem::float0;
  int ivfix = fem::int0;
  int nfron = fem::int0;
  int kfron = fem::int0;
  int kount = fem::int0;
  int nloca = fem::int0;
  int idofn = fem::int0;
  int kvfix = fem::int0;
  int jvfix = fem::int0;
  int keror = fem::int0;
  static const char* format_905 =
    "(/,/,' MAXIMUM FRONTWIDTH ENCOUNTERED =',i5)";
  static const char* format_910 =
    "(/,/,' *** DIAGNOSIS BY CHECK2, ERROR',i3,6x,' ASSOCIATED NUMBER',i5)";
  //C*********************************************************************
  //C
  //C*** THIS SUBROUTINE CHECKS THE REMAINDER OF THE INPUT DATA
  //C
  //C*********************************************************************
  //C
  //C*** CHECK AGAINST TWO IDENTICAL NOONZERO NODAL COORDINATES
  //C
  FEM_DO_SAFE(ieror, 13, 24) {
    neror(ieror) = 0;
  }
  FEM_DO_SAFE(ielem, 1, nelem) {
    ndfro(ielem) = 0;
  }
  FEM_DO_SAFE(ipoin, 2, npoin) {
    kpoin = ipoin - 1;
    FEM_DO_SAFE(jpoin, 1, kpoin) {
      FEM_DO_SAFE(idime, 1, 2) {
        if (coord(ipoin, idime) != coord(jpoin, idime)) {
          goto statement_30;
        }
      }
      neror(13)++;
      statement_30:;
    }
    //C
    //C*** CHECK THE LIST OF ELEMENT PROPERTY NUMBERS
    //C
    FEM_DO_SAFE(ielem, 1, nelem) {
      if (matno(ielem) <= 0 || matno(ielem) > nmats) {
        neror(14)++;
      }
    }
  }
  //C
  //C*** CHECK FOR IMPOSSIBLE NODE NUMBERS
  //C
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(inode, 1, nnode) {
      if (lnods(ielem, inode) == 0) {
        neror(15)++;
      }
      if (lnods(ielem, inode) < 0 || lnods(ielem, inode) > npoin) {
        neror(16)++;
      }
    }
  }
  //C
  //C*** CHECK FOR ANY REPETITION OF A NODE NUMBER WITHIN AN ELEMENT
  //C
  FEM_DO_SAFE(ipoin, 1, npoin) {
    kstar = 0;
    FEM_DO_SAFE(ielem, 1, nelem) {
      kzero = 0;
      FEM_DO_SAFE(inode, 1, nnode) {
        if (lnods(ielem, inode) != ipoin) {
          goto statement_90;
        }
        kzero++;
        if (kzero > 1) {
          neror(17)++;
        }
        //C
        //C*** SEEK FIRST,LAST AND INTERMEDIATE APPEARANCES OF NODE IPOIN
        //C
        if (kstar != 0) {
          goto statement_80;
        }
        kstar = ielem;
        //C
        //C*** CALCULATE INCREASE OR DECREASE IN FRONTWIDTH AT EACH ELEMENT STAGE
        //C
        ndfro(ielem) += ndofn;
        statement_80:
        //C
        //C*** AND CHANGE THE SIGN OF THE LAST APPEARANCE OF EACH NODE
        //C
        klast = ielem;
        nlast = inode;
        statement_90:;
      }
    }
    if (kstar == 0) {
      goto statement_110;
    }
    if (klast < nelem) {
      ndfro(klast + 1) = ndfro(klast + 1) - ndofn;
    }
    lnods(klast, nlast) = -ipoin;
    goto statement_140;
    //C
    //C*** CHECK THAT COORDINATES FOR AN UNUSED NODE HAVE NOT BEEN SPECIFIED
    //C
    statement_110:
    write(6, "(/,' CHECK WHY NODE',i4,' NEVER APPEARS')"), ipoin;
    neror(18)++;
    sigma = 0.0f;
    FEM_DO_SAFE(idime, 1, 2) {
      sigma += fem::abs(coord(ipoin, idime));
    }
    if (sigma != 0.0f) {
      neror(19)++;
    }
    //C
    //C*** CHECK THAT AN UNUSED NODE NUMBER IS NOT A RESTRAINED NODE
    //C
    FEM_DO_SAFE(ivfix, 1, nvfix) {
      if (nofix(ivfix) == ipoin) {
        neror(20)++;
      }
    }
    statement_140:;
  }
  //C
  //C*** CALCULATE THE LARGEST FRONTWIDTH
  //C
  nfron = 0;
  kfron = 0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    nfron += ndfro(ielem);
    if (nfron > kfron) {
      kfron = nfron;
    }
  }
  write(6, format_905), kfron;
  //write(6, format_905), kfron;
  if (kfron > mfron) {
    neror(21) = 1;
  }
  //C
  //C*** CONTINUE CHECKING THE DATA FOR THE FIXED VALUES
  //C
  FEM_DO_SAFE(ivfix, 1, nvfix) {
    if (nofix(ivfix) <= 0 || nofix(ivfix) > npoin) {
      neror(22)++;
    }
    kount = 0;
    nloca = (nofix(ivfix) - 1) * ndofn;
    FEM_DO_SAFE(idofn, 1, ndofn) {
      nloca++;
      if (iffix(nloca) > 0) {
        kount = 1;
      }
    }
    if (kount == 0) {
      neror(23)++;
    }
    kvfix = ivfix - 1;
    FEM_DO_SAFE(jvfix, 1, kvfix) {
      if (ivfix != 1 && nofix(ivfix) == nofix(jvfix)) {
        neror(24)++;
      }
    }
  }
  keror = 0;
  FEM_DO_SAFE(ieror, 13, 24) {
    if (neror(ieror) == 0) {
      goto statement_180;
    }
    keror = 1;
    write(6, format_910), ieror, neror(ieror);
    //write(6, format_910), ieror, neror(ieror);
    statement_180:;
  }
  if (keror != 0) {
    goto statement_200;
  }
  //C
  //C*** RETURN ALL NODAL CONNECTION NUMBERS TO POSITIVE VALUES
  //C
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(inode, 1, nnode) {
      lnods(ielem, inode) = fem::iabs(lnods(ielem, inode));
    }
  }
  return;
  statement_200:
  echo(cmn);
}

//C
//C(Page 17
void
nodexy(
  arr_ref<float, 2> coord,
  arr_cref<int, 2> lnods,
  int const& melem,
  int const& mpoin,
  int const& nelem,
  int const& nnode)
{
  coord(dimension(mpoin, 2));
  lnods(dimension(melem, 9));
  int ielem = fem::int0;
  int nnod1 = fem::int0;
  int inode = fem::int0;
  int nodst = fem::int0;
  int igash = fem::int0;
  int nodfn = fem::int0;
  int midpt = fem::int0;
  int nodmd = fem::int0;
  float total = fem::float0;
  int kount = fem::int0;
  int lnode = fem::int0;
  int lnod1 = fem::int0;
  int lnod3 = fem::int0;
  int lnod5 = fem::int0;
  int lnod7 = fem::int0;
  //C**********************************************************************
  //C
  //C****  THIS SUBROUTINE INTERPOLATES THE MID SIDE NODES OF STRAIGHT
  //C       SIDES OF ELEMENTS AND THE CENTRAL NODE OF 9 NODED ELEMENTS
  //C
  //C**********************************************************************
  if (nnode == 4) {
    return;
  }
  //C
  //C*** LOOP OVER EACH ELEMENT
  //C
  FEM_DO_SAFE(ielem, 1, nelem) {
    //C
    //C*** LOOP OVER EACH ELEMENT EDGE
    //C
    nnod1 = 9;
    if (nnode == 8) {
      nnod1 = 7;
    }
    FEM_DOSTEP(inode, 1, nnod1, 2) {
      if (inode == 9) {
        goto statement_50;
      }
      //C
      //C*** COMPUTE THE NODE NUMBER OF THE FIRST NODE
      //C
      nodst = lnods(ielem, inode);
      igash = inode + 2;
      if (igash > 8) {
        igash = 1;
      }
      //C
      //C*** COMPUTE THE NODE NUMBER OF THE LAST NODE
      //C
      nodfn = lnods(ielem, igash);
      midpt = inode + 1;
      //C
      //C*** COMPUTE THE NODE NUMBER OF THE INTERMEDIATE NODE
      //C
      nodmd = lnods(ielem, midpt);
      total = fem::abs(coord(nodmd, 1)) + fem::abs(coord(nodmd, 2));
      //C
      //C*** IF THE COORDINATES OF THE INTERMEDIATE NODE ARE BOTH
      //C    ZERO INTERPOLATE BY A STRAIGHT LINE
      //C
      if (total > 0.0f) {
        goto statement_20;
      }
      kount = 1;
      statement_10:
      coord(nodmd, kount) = (coord(nodst, kount) + coord(nodfn, kount)) / 2.0f;
      kount++;
      if (kount == 2) {
        goto statement_10;
      }
      statement_20:;
    }
    goto statement_30;
    statement_50:
    lnode = lnods(ielem, inode);
    total = fem::abs(coord(lnode, 1)) + fem::abs(coord(lnode, 2));
    if (total > 0.0f) {
      goto statement_30;
    }
    lnod1 = lnods(ielem, 1);
    lnod3 = lnods(ielem, 3);
    lnod5 = lnods(ielem, 5);
    lnod7 = lnods(ielem, 7);
    kount = 1;
    statement_40:
    coord(lnode, kount) = (coord(lnod1, kount) + coord(lnod3,
      kount) + coord(lnod5, kount) + coord(lnod7, kount)) / 4.0f;
    kount++;
    if (kount == 2) {
      goto statement_40;
    }
    statement_30:;
  }
}

//C
//C(Page 179)
void
gaussq(
  int const& ngaus,
  arr_ref<float> posgp,
  arr_ref<float> weigp)
{
  posgp(dimension(4));
  weigp(dimension(4));
  int kgaus = fem::int0;
  int igash = fem::int0;
  int jgash = fem::int0;
  //C*******************************************************************
  //C
  //C *** SETS UP THE GAUSS-LEGENDRE INTEGRATION CONSTANTS
  //C
  //C*******************************************************************
  if (ngaus > 2) {
    goto statement_4;
  }
  posgp(1) = -0.577350269189626f;
  weigp(1) = 1.0f;
  goto statement_6;
  statement_4:
  posgp(1) = -0.774596669241483f;
  posgp(2) = 0.0f;
  weigp(1) = 0.5555555555555556f;
  weigp(2) = 0.8888888888888889f;
  statement_6:
  kgaus = ngaus / 2;
  FEM_DO_SAFE(igash, 1, kgaus) {
    jgash = ngaus + 1 - igash;
    posgp(jgash) = -posgp(igash);
    weigp(jgash) = weigp(igash);
  }
}

//C
//C(Page 180)
void
sfr2(
  arr_ref<float, 2> deriv,
  float const& etasp,
  float const& exisp,
  int const& nnode,
  arr_ref<float> shape)
{
  deriv(dimension(2, 9));
  shape(dimension(9));
  float s = fem::float0;
  float t = fem::float0;
  float st = fem::float0;
  float s2 = fem::float0;
  float t2 = fem::float0;
  float ss = fem::float0;
  float tt = fem::float0;
  float sst = fem::float0;
  float stt = fem::float0;
  float st2 = fem::float0;
  float s1 = fem::float0;
  float t1 = fem::float0;
  float s9 = fem::float0;
  float t9 = fem::float0;
  //C******************************************************************
  //C
  //C *** THIS SUBROUTINUE EVALUATES SHAPE FUNCTIONS AND THEIR
  //C     DERIVERTIVES FOR LINEAR QUADRATIC LAGRANGIAN AND
  //C     SERENDIPITY ISOPARAMETRIC 2-D ELEMENTS
  //C
  //C*******************************************************************
  s = exisp;
  t = etasp;
  if (nnode > 4) {
    goto statement_10;
  }
  st = s * t;
  //C
  //C *** SHAPE FUNCTION FOR 4 NODED ELEMENT
  //C
  shape(1) = (1 - t - s + st) * 0.25f;
  shape(2) = (1 - t + s - st) * 0.25f;
  shape(3) = (1 + t + s + st) * 0.25f;
  shape(4) = (1 + t - s - st) * 0.25f;
  //C
  //C *** SHAPE FUNCTION DERIVERTIVES
  //C
  deriv(1, 1) = (-1 + t) * 0.25f;
  deriv(1, 2) = (+1 - t) * 0.25f;
  deriv(1, 3) = (+1 + t) * 0.25f;
  deriv(1, 4) = (-1 - t) * 0.25f;
  deriv(2, 1) = (-1 + s) * 0.25f;
  deriv(2, 2) = (-1 - s) * 0.25f;
  deriv(2, 3) = (+1 + s) * 0.25f;
  deriv(2, 4) = (+1 - s) * 0.25f;
  return;
  statement_10:
  if (nnode > 8) {
    goto statement_30;
  }
  s2 = s * 2.0f;
  t2 = t * 2.0f;
  ss = s * s;
  tt = t * t;
  st = s * t;
  sst = s * s * t;
  stt = s * t * t;
  st2 = s * t * 2.0f;
  //C
  //C *** SHAPE FUNCTIONS FOR 8 NODED ELEMENT
  //C
  shape(1) = (-1.0f + st + ss + tt - sst - stt) / 4.0f;
  shape(2) = (1.0f - t - ss + sst) / 2.0f;
  shape(3) = (-1.0f - st + ss + tt - sst + stt) / 4.0f;
  shape(4) = (1.0f + s - tt - stt) / 2.0f;
  shape(5) = (-1.0f + st + ss + tt + sst + stt) / 4.0f;
  shape(6) = (1.0f + t - ss - sst) / 2.0f;
  shape(7) = (-1.0f - st + ss + tt + sst - stt) / 4.0f;
  shape(8) = (1.0f - s - tt + stt) / 2.0f;
  //C
  //C *** SHAPE FUNCTION DERIVERTIVES
  //C
  deriv(1, 1) = (t + s2 - st2 - tt) / 4.0f;
  deriv(1, 2) = -s + st;
  deriv(1, 3) = (-t + s2 - st2 + tt) / 4.0f;
  deriv(1, 4) = (1.0f - tt) / 2.0f;
  deriv(1, 5) = (t + s2 + st2 + tt) / 4.0f;
  deriv(1, 6) = -s - st;
  deriv(1, 7) = (-t + s2 + st2 - tt) / 4.0f;
  deriv(1, 8) = (-1.0f + tt) / 2.0f;
  deriv(2, 1) = (s + t2 - ss - st2) / 4.0f;
  deriv(2, 2) = (-1.0f + ss) / 2.0f;
  deriv(2, 3) = (-s + t2 - ss + st2) / 4.0f;
  deriv(2, 4) = -t - st;
  deriv(2, 5) = (s + t2 + ss + st2) / 4.0f;
  deriv(2, 6) = (1.0f - ss) / 2.0f;
  deriv(2, 7) = (-s + t2 + ss - st2) / 4.0f;
  deriv(2, 8) = -t + st;
  return;
  statement_30:
  ss = s * s;
  st = s * t;
  tt = t * t;
  s1 = s + 1.0f;
  t1 = t + 1.0f;
  s2 = s * 2.0f;
  t2 = t * 2.0f;
  s9 = s - 1.0f;
  t9 = t - 1.0f;
  //C
  //C *** SHAPE FUNCTION FOR 9 NODED ELEMENT
  //C
  shape(1) = 0.25f * s9 * st * t9;
  shape(2) = 0.5f * (1.0f - ss) * t * t9;
  shape(3) = 0.25f * s1 * st * t9;
  shape(4) = 0.5f * s * s1 * (1.0f - tt);
  shape(5) = 0.25f * s1 * st * t1;
  shape(6) = 0.5f * (1.0f - ss) * t * t1;
  shape(7) = 0.25f * s9 * st * t1;
  shape(8) = 0.5f * s * s9 * (1.0f - tt);
  shape(9) = (1.0f - ss) * (1.0f - tt);
  //C
  //C *** SHAPR FUNCTION DERIVERTIVES
  //C
  deriv(1, 1) = 0.25f * t * t9 * (-1.0f + s2);
  deriv(1, 2) = -st * t9;
  deriv(1, 3) = 0.25f * (1.0f + s2) * t * t9;
  deriv(1, 4) = 0.5f * (1.0f + s2) * (1.0f - tt);
  deriv(1, 5) = 0.25f * (1.0f + s2) * t * t1;
  deriv(1, 6) = -st * t1;
  deriv(1, 7) = 0.25f * (-1.0f + s2) * t * t1;
  deriv(1, 8) = 0.5f * (-1.0f + s2) * (1.0f - tt);
  deriv(1, 9) = -s2 * (1.0f - tt);
  deriv(2, 1) = 0.25f * (-1.0f + t2) * s * s9;
  deriv(2, 2) = 0.5f * (1.0f - ss) * (-1.0f + t2);
  deriv(2, 3) = 0.25f * s * s1 * (-1.0f + t2);
  deriv(2, 4) = -st * s1;
  deriv(2, 5) = 0.25f * s * s1 * (1.0f + t2);
  deriv(2, 6) = 0.5f * (1.0f - ss) * (1.0f + t2);
  deriv(2, 7) = 0.25f * s * s9 * (1.0f + t2);
  deriv(2, 8) = -st * s9;
  deriv(2, 9) = -t2 * (1.0f - ss);
}

//C
//C(Page 182
void
jacob2(
  common& cmn,
  arr_ref<float, 2> cartd,
  arr_cref<float, 2> deriv,
  float& djacb,
  arr_cref<float, 2> elcod,
  arr_ref<float, 2> gpcod,
  int const& ielem,
  int const& kgasp,
  int const& nnode,
  arr_cref<float> shape)
{
  cartd(dimension(2, 9));
  deriv(dimension(2, 9));
  elcod(dimension(2, 9));
  gpcod(dimension(2, 9));
  shape(dimension(9));
  common_write write(cmn);
  int idime = fem::int0;
  int inode = fem::int0;
  int jdime = fem::int0;
  arr_2d<2, 2, float> xjacm(fem::fill0);
  arr_2d<2, 2, float> xjaci(fem::fill0);
  //C********************************************************************
  //C
  //C *** THIS SUBROUTINE EVALUATES THE JACOBIAN MATRIX AND THE
  //C     CARTESIAN SHAPE FUNCTION DERIVERTIVES
  //C
  //C********************************************************************
  //C
  //C *** CALCULATE COORDINATES OF SAMPLING POINT
  //C
  FEM_DO_SAFE(idime, 1, 2) {
    gpcod(idime, kgasp) = 0.0f;
    FEM_DO_SAFE(inode, 1, nnode) {
      gpcod(idime, kgasp) += elcod(idime, inode) * shape(inode);
    }
  }
  //C
  //C *** CREATE JACOBIAN MATRIX XJACM
  //C
  FEM_DO_SAFE(idime, 1, 2) {
    FEM_DO_SAFE(jdime, 1, 2) {
      xjacm(idime, jdime) = 0.0f;
      FEM_DO_SAFE(inode, 1, nnode) {
        xjacm(idime, jdime) += deriv(idime, inode) * elcod(jdime, inode);
      }
    }
  }
  //C
  //C *** CALCULATE DETERMINANT AND INVERSE OF JACOBIAN MATRIX
  //C
  djacb = xjacm(1, 1) * xjacm(2, 2) - xjacm(1, 2) * xjacm(2, 1);
  switch (fem::if_arithmetic(djacb)) {
    case -1: goto statement_6;
    case  0: goto statement_6;
    default: goto statement_8;
  }
  statement_6:
  write(6,
    "(/,/,' PROGRAM HALTED IN SUBROUTINE JACOB2',/,11x,"
    "' ZERO OR NEGATIVE AREA',/,10x,' ELEMENT NUMBER',i5)"),
    ielem;
  FEM_STOP(0);
  statement_8:
  xjaci(1, 1) = xjacm(2, 2) / djacb;
  xjaci(2, 2) = xjacm(1, 1) / djacb;
  xjaci(1, 2) = -xjacm(1, 2) / djacb;
  xjaci(2, 1) = -xjacm(2, 1) / djacb;
  //C
  //C *** CALCULATE CARTESIAN DERIVERTIVES
  //C
  FEM_DO_SAFE(idime, 1, 2) {
    FEM_DO_SAFE(inode, 1, nnode) {
      cartd(idime, inode) = 0.0f;
      FEM_DO_SAFE(jdime, 1, 2) {
        cartd(idime, inode) += xjaci(idime, jdime) * deriv(jdime, inode);
      }
    }
  }
}

//C
//C(P
void
bmatps(
  arr_ref<float, 2> bmatx,
  arr_cref<float, 2> cartd,
  int const& nnode,
  arr_cref<float> shape,
  arr_cref<float, 2> gpcod,
  int const& ntype,
  int const& kgasp)
{
  bmatx(dimension(4, 18));
  cartd(dimension(2, 9));
  shape(dimension(9));
  gpcod(dimension(2, 9));
  int ngash = fem::int0;
  int inode = fem::int0;
  int mgash = fem::int0;
  //C**********************************************************************
  //C
  //C *** THIS SUBROUTINE EVALUATES THE STRAIN-DISPLACEMENT MATRIX
  //C
  //C**********************************************************************
  ngash = 0;
  FEM_DO_SAFE(inode, 1, nnode) {
    mgash = ngash + 1;
    ngash = mgash + 1;
    bmatx(1, mgash) = cartd(1, inode);
    bmatx(1, ngash) = 0.0f;
    bmatx(2, mgash) = 0.0f;
    bmatx(2, ngash) = cartd(2, inode);
    bmatx(3, mgash) = cartd(2, inode);
    bmatx(3, ngash) = cartd(1, inode);
    if (ntype != 3) {
      goto statement_10;
    }
    bmatx(4, mgash) = shape(inode) / gpcod(1, kgasp);
    bmatx(4, ngash) = 0.0f;
    statement_10:;
  }
}

//C
//C(Page 192)
void
modps(
  arr_ref<float, 2> dmatx,
  int const& lprop,
  int const& mmats,
  int const& ntype,
  arr_cref<float, 2> props)
{
  dmatx(dimension(4, 4));
  props(dimension(mmats, 7));
  float young = fem::float0;
  float poiss = fem::float0;
  int istr1 = fem::int0;
  int jstr1 = fem::int0;
  float identifier_const = fem::float0;
  float conss = fem::float0;
  //C*********************************************************************
  //C
  //C *** THIS SUBROUTINE EVALUATES THE D-MATRIX
  //C
  //C*********************************************************************
  young = props(lprop, 1);
  poiss = props(lprop, 2);
  FEM_DO_SAFE(istr1, 1, 4) {
    FEM_DO_SAFE(jstr1, 1, 4) {
      dmatx(istr1, jstr1) = 0.0f;
    }
  }
  if (ntype != 1) {
    goto statement_4;
  }
  //C
  //C *** D MATRIX FOR PLANE STRESS CASE
  //C
  identifier_const = young / (1.0f - poiss * poiss);
  dmatx(1, 1) = identifier_const;
  dmatx(2, 2) = identifier_const;
  dmatx(1, 2) = identifier_const * poiss;
  dmatx(2, 1) = identifier_const * poiss;
  dmatx(3, 3) = (1.0f - poiss) * identifier_const / 2.0f;
  return;
  statement_4:
  if (ntype != 2) {
    goto statement_6;
  }
  //C
  //C *** D MATRIX FOR PLANE STRAIN CASE
  //C
  identifier_const = young * (1.0f - poiss) / ((1.0f + poiss) * (
    1.0f - 2.0f * poiss));
  dmatx(1, 1) = identifier_const;
  dmatx(2, 2) = identifier_const;
  dmatx(1, 2) = identifier_const * poiss / (1.0f - poiss);
  dmatx(2, 1) = identifier_const * poiss / (1.0f - poiss);
  dmatx(3, 3) = (1.0f - 2.0f * poiss) * identifier_const / (2.0f * (
    1.0f - poiss));
  return;
  statement_6:
  if (ntype != 3) {
    goto statement_8;
  }
  //C
  //C *** D MATRIX FOR AXISYMMETRIC CASE
  //C
  identifier_const = young * (1.0f - poiss) / ((1.0f + poiss) * (
    1.0f - 2.0f * poiss));
  conss = poiss / (1.0f - poiss);
  dmatx(1, 1) = identifier_const;
  dmatx(2, 2) = identifier_const;
  dmatx(3, 3) = identifier_const * (1.0f - 2.0f * poiss) / (2.0f * (
    1.0f - poiss));
  dmatx(1, 2) = identifier_const * conss;
  dmatx(1, 4) = identifier_const * conss;
  dmatx(2, 1) = identifier_const * conss;
  dmatx(2, 4) = identifier_const * conss;
  dmatx(4, 1) = identifier_const * conss;
  dmatx(4, 2) = identifier_const * conss;
  dmatx(4, 4) = identifier_const;
  statement_8:;
}

//C
//C(Page
void
dbe(
  arr_cref<float, 2> bmatx,
  arr_ref<float, 2> dbmat,
  arr_cref<float, 2> dmatx,
  int const& mevab,
  int const& nevab,
  int const& nstre,
  int const& nstr1)
{
  bmatx(dimension(nstr1, mevab));
  dbmat(dimension(nstr1, mevab));
  dmatx(dimension(nstr1, nstr1));
  //C**********************************************************************
  //C
  //C *** THIS SUBROUTINE MULTIPLIES THE D-MATRIX BY B-MATRIX
  //C
  //C**********************************************************************
  int istre = fem::int0;
  int ievab = fem::int0;
  int jstre = fem::int0;
  FEM_DO_SAFE(istre, 1, nstre) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      dbmat(istre, ievab) = 0.0f;
      FEM_DO_SAFE(jstre, 1, nstre) {
        dbmat(istre, ievab) += dmatx(istre, jstre) * bmatx(jstre, ievab);
      }
    }
  }
}

//C$DEBUG
//C$LARGE
//C(Page 213)
void
conver(
  common& cmn,
  arr_ref<float, 2> eload,
  int const& iiter,
  arr_cref<int, 2> lnods,
  int const& melem,
  int const& mevab,
  int const& mtotv,
  int& nchek,
  int const& ndofn,
  int const& nelem,
  int const& nevab,
  int const& nnode,
  int const& ntotv,
  float& pvalu,
  arr_ref<float> stfor,
  arr_cref<float, 2> tload,
  arr_ref<float> tofor,
  float const& toler)
{
  eload(dimension(melem, mevab));
  lnods(dimension(melem, 9));
  stfor(dimension(mtotv));
  tload(dimension(melem, mevab));
  tofor(dimension(mtotv));
  common_write write(cmn);
  float resid = fem::float0;
  float retot = fem::float0;
  float remax = fem::float0;
  int itotv = fem::int0;
  int ielem = fem::int0;
  int kevab = fem::int0;
  int inode = fem::int0;
  int locno = fem::int0;
  int idofn = fem::int0;
  int nposi = fem::int0;
  float refor = fem::float0;
  float agash = fem::float0;
  int ievab = fem::int0;
  float ratio = fem::float0;
  static const char* format_30 =
    "('0',3x,'CONVERGENCE CODE =',i4,3x,'NORM OF RESIDUAL SUMRATIO =,',e14.6,"
    "3x,'MAXIMUM RESIDUAL =',e14.6)";
  //C**********************************************************************
  //C
  //C*** THIS SUBROUTINE CHECKS FOR CONVERGENCE OF THE ITERATION PROCESS
  //C
  //C**********************************************************************
  nchek = 0;
  resid = 0.0f;
  retot = 0.0f;
  remax = 0.0f;
  FEM_DO_SAFE(itotv, 1, ntotv) {
    stfor(itotv) = 0.0f;
    tofor(itotv) = 0.0f;
  }
  FEM_DO_SAFE(ielem, 1, nelem) {
    kevab = 0;
    FEM_DO_SAFE(inode, 1, nnode) {
      locno = fem::iabs(lnods(ielem, inode));
      FEM_DO_SAFE(idofn, 1, ndofn) {
        kevab++;
        nposi = (locno - 1) * ndofn + idofn;
        stfor(nposi) += eload(ielem, kevab);
        tofor(nposi) += tload(ielem, kevab);
      }
    }
  }
  FEM_DO_SAFE(itotv, 1, ntotv) {
    refor = tofor(itotv) - stfor(itotv);
    resid += refor * refor;
    retot += tofor(itotv) * tofor(itotv);
    agash = fem::abs(refor);
    if (agash > remax) {
      remax = agash;
    }
  }
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      eload(ielem, ievab) = tload(ielem, ievab) - eload(ielem, ievab);
    }
  }
  resid = fem::sqrt(resid);
  retot = fem::sqrt(retot);
  ratio = 100.0f * resid / retot;
  if (ratio > toler) {
    nchek = 1;
  }
  if (iiter == 1) {
    goto statement_20;
  }
  if (ratio > pvalu) {
    nchek = 999;
  }
  statement_20:
  pvalu = ratio;
  write(6, format_30), nchek, ratio, remax;
  //write(6, format_30), nchek, ratio, remax;
}

//C$DEBUG
//C$LARGE
//C(Page 238)
void
dimen(
  int& mbufa,
  int& melem,
  int& mevab,
  int& mfron,
  int& mmats,
  int& mpoin,
  int& mstif,
  int& mtotg,
  int& mtotv,
  int& mvfix,
  int& ndofn,
  int& nprop,
  int const& /* nstre */)
{
  //C**********************************************************************
  //C
  //C*** THIS SUBROUTINE PRESETS VARIABLES ASSOCIATED WITH DYNAMIC
  //C    DIMENSIONING
  //C
  //C**********************************************************************
  mbufa = 10;
  melem = 40;
  mfron = 80;
  mmats = 5;
  mpoin = 150;
  mstif = (mfron * mfron - mfron) / 2.0f + mfron;
  mtotg = melem * 9;
  ndofn = 2;
  mtotv = mpoin * ndofn;
  mvfix = 30;
  nprop = 7;
  mevab = ndofn * 9;
}

//C$DEBUG
//C$LARGE
//C(Page 243)
void
flowpl(
  arr_cref<float> avect,
  float& abeta,
  arr_ref<float> dvect,
  int const& ntype,
  arr_cref<float, 2> props,
  int const& lprop,
  int const& nstr1,
  int const& mmats)
{
  avect(dimension(4));
  dvect(dimension(4));
  props(dimension(mmats, 7));
  float young = fem::float0;
  float poiss = fem::float0;
  float hards = fem::float0;
  float fmul1 = fem::float0;
  float fmul2 = fem::float0;
  float fmul3 = fem::float0;
  float denom = fem::float0;
  int istr1 = fem::int0;
  //C******************************************************************
  //C
  //C****THIS SUBROUTINE EVALUATES THE PLASTIC D VECTOR
  //C
  //C******************************************************************
  young = props(lprop, 1);
  poiss = props(lprop, 2);
  hards = props(lprop, 6);
  fmul1 = young / (1.0f + poiss);
  if (ntype == 1) {
    goto statement_60;
  }
  fmul2 = young * poiss * (avect(1) + avect(2) + avect(4)) / ((1.0f +
    poiss) * (1.0f - 2.0f * poiss));
  dvect(1) = fmul1 * avect(1) + fmul2;
  dvect(2) = fmul1 * avect(2) + fmul2;
  dvect(3) = 0.5f * avect(3) * young / (1.0f + poiss);
  dvect(4) = fmul1 * avect(4) + fmul2;
  goto statement_70;
  statement_60:
  fmul3 = young * poiss * (avect(1) + avect(2)) / (1.0f - poiss * poiss);
  dvect(1) = fmul1 * avect(1) + fmul3;
  dvect(2) = fmul1 * avect(2) + fmul3;
  dvect(3) = 0.5f * avect(3) * young / (1.0f + poiss);
  dvect(4) = fmul1 * avect(4) + fmul3;
  statement_70:
  denom = hards;
  FEM_DO_SAFE(istr1, 1, nstr1) {
    denom += avect(istr1) * dvect(istr1);
  }
  abeta = 1.0f / denom;
}

int nfunc(int i,int j) {
	return  i + (j * j - j) / 2;
}
//C$DEBUG
//C$LARGE
//C(Page 194)
void
front(
  common& cmn,
  arr_ref<float> asdis,
  arr_cref<float, 2> eload,
  arr_ref<float> eqrhs,
  arr_ref<float, 2> equat,
  arr_ref<float, 2> estif,
  arr_ref<float> fixed,
  arr_cref<int> iffix,
  int const& iincs,
  int const& iiter,
  arr_ref<float> gload,
  arr_ref<float> gstif,
  arr_ref<int> locel,
  arr_ref<int, 2> lnods,
  int const& kresl,
  int const& mbufa,
  int const& melem,
  int const& mevab,
  int const& mfron,
  int const& mstif,
  int const& mtotv,
  int const& mvfix,
  arr_ref<int> nacva,
  arr_ref<int> namev,
  arr_ref<int> ndest,
  int const& ndofn,
  int const& nelem,
  int const& nevab,
  int const& nnode,
  arr_cref<int> /* nofix */,
  arr_ref<int> npivo,
  int const& npoin,
  int const& ntotv,
  arr_ref<float> tdisp,
  arr_ref<float, 2> tload,
  arr_ref<float, 2> treac,
  arr_ref<float> vecrv)
{
  asdis(dimension(mtotv));
  eload(dimension(melem, mevab));
  eqrhs(dimension(mbufa));
  equat(dimension(mfron, mbufa));
  estif(dimension(mevab, mevab));
  fixed(dimension(mtotv));
  iffix(dimension(mtotv));
  gload(dimension(mfron));
  gstif(dimension(mstif));
  locel(dimension(mevab));
  lnods(dimension(melem, 9));
  nacva(dimension(mfron));
  namev(dimension(mbufa));
  ndest(dimension(mevab));
  npivo(dimension(mbufa));
  tdisp(dimension(mtotv));
  tload(dimension(melem, mevab));
  treac(dimension(mvfix, ndofn));
  vecrv(dimension(mfron));
  common_read read(cmn);
  common_write write(cmn);
  int i = fem::int0;
  int j = fem::int0;
  int ipoin = fem::int0;
  int klast = fem::int0;
  int ielem = fem::int0;
  int inode = fem::int0;
  int nlast = fem::int0;
  int ibufa = fem::int0;
  int istif = fem::int0;
  int ifron = fem::int0;
  int nbufa = fem::int0;
  int nfron = fem::int0;
  int kelva = fem::int0;
  int kevab = fem::int0;
  int idofn = fem::int0;
  int nposi = fem::int0;
  int locno = fem::int0;
  int ievab = fem::int0;
  int nikno = fem::int0;
  int kexis = fem::int0;
  int idest = fem::int0;
  int jevab = fem::int0;
  int jdest = fem::int0;
  int ngash = fem::int0;
  int ngish = fem::int0;
  int jfron = fem::int0;
  int nloca = fem::int0;
  float pivot = fem::float0;
  float cureq = fem::float0;
  int lfron = fem::int0;
  int ielva = fem::int0;
  int itotv = fem::int0;
  int kboun = fem::int0;
  int ngush = fem::int0;
  int mgash = fem::int0;
  //C**********************************************************************
  //C                                                                     *
  //C*** THIS SUBROUTINE UNDERTAKES EQUATION SOLUTION BY THE FRONTAL      *
  //C    METHOD                                                           *
  //C                                                                                                    *
  //C**********************************************************************
  //C
  //C**CHANGE THE SIGN OF THE LAST APPEARANCE OF EACH NODE **
  //C
  if (iincs > 1 || iiter > 1) {
    goto statement_455;
  }
  FEM_DO_SAFE(ipoin, 1, npoin) {
    klast = 0;
    FEM_DO_SAFE(ielem, 1, nelem) {
      FEM_DO_SAFE(inode, 1, nnode) {
        if (lnods(ielem, inode) != ipoin) {
          goto statement_120;
        }
        klast = ielem;
        nlast = inode;
        statement_120:;
      }
    }
    if (klast != 0) {
      lnods(klast, nlast) = -ipoin;
    }
  }
  statement_455:
  //C
  //C         ** START BY INITIALIZING EVERYTHING THAT MATTERS TO ZERO **
  //C
  FEM_DO_SAFE(ibufa, 1, mbufa) {
    eqrhs(ibufa) = 0.0f;
  }
  FEM_DO_SAFE(istif, 1, mstif) {
    gstif(istif) = 0.0f;
  }
  FEM_DO_SAFE(ifron, 1, mfron) {
    gload(ifron) = 0.0f;
    vecrv(ifron) = 0.0f;
    nacva(ifron) = 0;
    FEM_DO_SAFE(ibufa, 1, mbufa) {
      equat(ifron, ibufa) = 0.0f;
    }
  }
  //C
  //C        ** AND PREPARE FOR DISC READING AND WRITING OPERATIONS **
  //C
  nbufa = 0;
  if (kresl > 1) {
    nbufa = mbufa;
  }
  cmn.io.rewind(1);
  cmn.io.rewind(2);
  cmn.io.rewind(3);
  cmn.io.rewind(4);
  cmn.io.rewind(8);
  //C
  //C        ** ENTER MAIN ELEMENT ASSEMBLY-REDUCTION LOOP ***
  //C
  nfron = 0;
  kelva = 0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    if (kresl > 1) {
      goto statement_400;
    }
    kevab = 0;
    read(1, fem::unformatted), estif;
    FEM_DO_SAFE(inode, 1, nnode) {
      FEM_DO_SAFE(idofn, 1, ndofn) {
        nposi = (inode - 1) * ndofn + idofn;
        locno = lnods(ielem, inode);
        if (locno > 0) {
          locel(nposi) = (locno - 1) * ndofn + idofn;
        }
        if (locno < 0) {
          locel(nposi) = (locno + 1) * ndofn - idofn;
        }
      }
    }
    //C
    //C        ** START BY LOOKING FOR EXISTING DESTINATIONS *******
    //C
    FEM_DO_SAFE(ievab, 1, nevab) {
      nikno = fem::iabs(locel(ievab));
      kexis = 0;
      FEM_DO_SAFE(ifron, 1, nfron) {
        if (nikno != nacva(ifron)) {
          goto statement_180;
        }
        kevab++;
        kexis = 1;
        ndest(kevab) = ifron;
        statement_180:;
      }
      if (kexis != 0) {
        goto statement_210;
      }
      //C
      //C***  WE NOW SEEK NEW EMPTY PLACES FOR DESTINATION VECTOR ***
      //C
      FEM_DO_SAFE(ifron, 1, mfron) {
        if (nacva(ifron) != 0) {
          goto statement_190;
        }
        nacva(ifron) = nikno;
        kevab++;
        ndest(kevab) = ifron;
        goto statement_200;
        statement_190:;
      }
      //C
      //C         THE NEW PLACES MAY DEMAND AN INCREASE IN CURRENT FRONTWIDTH**
      //C
      statement_200:
      if (ndest(kevab) > nfron) {
        nfron = ndest(kevab);
      }
      statement_210:;
    }
    write(8, fem::unformatted), locel, ndest, nacva, nfron;
    statement_400:
    if (kresl > 1) {
      read(8, fem::unformatted), locel, ndest, nacva, nfron;
    }
    //C
    //C***  ASSEMBLE ELEMENT LOADS
    //C
    FEM_DO_SAFE(ievab, 1, nevab) {
      idest = ndest(ievab);
      gload(idest) += eload(ielem, ievab);
      //C
      //C***  ASSEMBLE THE ELEMENT STIFFNESSES BUT NOT IN RESOLUTION***
      //C
      if (kresl > 1) {
        goto statement_402;
      }
      FEM_DO_SAFE(jevab, 1, ievab) {
        jdest = ndest(jevab);
        ngash = nfunc(idest, jdest);
        ngish = nfunc(jdest, idest);
        if (jdest >= idest) {
          gstif(ngash) += estif(ievab, jevab);
        }
        if (jdest < idest) {
          gstif(ngish) += estif(ievab, jevab);
        }
      }
      statement_402:;
    }
    //C
    //C*** RE-EXAMINE EACH ELEMENT NODE,TO ENQUIRE WHICH CAN BE ELIMINATED
    //C
    FEM_DO_SAFE(ievab, 1, nevab) {
      nikno = -locel(ievab);
      if (nikno <= 0) {
        goto statement_310;
      }
      //C
      //C***  FIND POSITIONS OF VARIABLES READY FOR ELIMINATION ****
      //C
      FEM_DO_SAFE(ifron, 1, nfron) {
        if (nacva(ifron) != nikno) {
          goto statement_300;
        }
        nbufa++;
        //C
        //C*** WRITE EQUITIONS TO DISC OR TO TAPE
        //C
        if (nbufa <= mbufa) {
          goto statement_406;
        }
        nbufa = 1;
        if (kresl > 1) {
          goto statement_408;
        }
		// zhangzhg
        write(2, fem::unformatted), equat, eqrhs, npivo, namev;
        goto statement_406;
        statement_408:
        write(4, fem::unformatted), eqrhs;
        read(2, fem::unformatted), equat, eqrhs, npivo, namev;
        statement_406:
        //C
        //C*** EXTRACT THE COEFFICIENTS OF THE NEW EQUATION FOR ELIMINATION
        //C
        if (kresl > 1) {
          goto statement_404;
        }
        FEM_DO_SAFE(jfron, 1, mfron) {
          if (ifron < jfron) {
            nloca = nfunc(ifron, jfron);
          }
          if (ifron >= jfron) {
            nloca = nfunc(jfron, ifron);
          }
          equat(jfron, nbufa) = gstif(nloca);
          gstif(nloca) = 0.0f;
        }
        statement_404:
        //C
        //C**** AND EXTRACT THE CORRESPONDING RIGHT HAND SIDES
        //C
        eqrhs(nbufa) = gload(ifron);
        gload(ifron) = 0.0f;
        kelva++;
        namev(nbufa) = nikno;
        npivo(nbufa) = ifron;
        //C
        //C***  DEAL WITH PIVOT
        //C
        pivot = equat(ifron, nbufa);
        if (pivot > 0.0f) {
          goto statement_235;
        }
        write(6,
          "('0',3x,'NEGATIVE OR ZERO PIVOT ENCOUNTERED FOR VARIABLE NO. ',i4,"
          "' OF VALUE ',e17.6)"),
          nikno, pivot;
        FEM_STOP(0);
        statement_235:
        equat(ifron, nbufa) = 0.0f;
        //C
        //C***  ENQUIRE WHETHER PRESENT VARIABLE IS FREE OR PRESCRIBED
        //C
        if (iffix(nikno) == 0) {
          goto statement_250;
        }
        //C
        //C*** DEAL WITH A PRESCRIBED DEFLECTION
        //C
        FEM_DO_SAFE(jfron, 1, nfron) {
          gload(jfron) = gload(jfron) - fixed(nikno) * equat(jfron, nbufa);
        }
        goto statement_280;
        //C
        //C*** ELIMINATE A FREE VARIABLE - DEAL WITH THE RIGHT HAND SIDE FIRST
        //C
        statement_250:
        FEM_DO_SAFE(jfron, 1, nfron) {
          gload(jfron) = gload(jfron) - equat(jfron, nbufa) * eqrhs(
            nbufa) / pivot;
          //C
          //C*** NOW DEAL WITH THE COEFFICIENTS IN CORE
          //C
          if (kresl > 1) {
            goto statement_418;
          }
          if (equat(jfron, nbufa) == 0.0f) {
            goto statement_270;
          }
          nloca = nfunc(0, jfron);
          cureq = equat(jfron, nbufa);
          FEM_DO_SAFE(lfron, 1, jfron) {
            ngash = lfron + nloca;
            gstif(ngash) = gstif(ngash) - cureq * equat(lfron, nbufa) / pivot;
          }
          statement_418:
          statement_270:;
        }
        statement_280:
        equat(ifron, nbufa) = pivot;
        //C
        //C*** RECORD THE NEW VACANT SPACE,AND REDUCE FRONTWIDTH IF POSSIBL
        //C
        nacva(ifron) = 0;
        goto statement_290;
        //C
        //C*** COMPLETE THE ELEMENT LOOP IN THE FORWARD ELEMINATION
        //C
        statement_300:;
      }
      statement_290:
      if (nacva(nfron) != 0) {
        goto statement_310;
      }
      nfron = nfron - 1;
      if (nfron > 0) {
        goto statement_290;
      }
      statement_310:;
    }
  }
  if (kresl == 1) {
    write(2, fem::unformatted), equat, eqrhs, npivo, namev;
  }
  cmn.io.backspace(2);
  //C
  //C*** ENTER BACK-SUBSROUTINE PHASE,LOOP BACKWORDS  THROUGH VARIABLES
  //C
  FEM_DO_SAFE(ielva, 1, kelva) {
    //C
    //C*** READ A NEW BLOCK OF EQUATIONS- IF NEEDED
    //C
    if (nbufa != 0) {
      goto statement_412;
    }
    cmn.io.backspace(2);
    read(2, fem::unformatted), equat, eqrhs, npivo, namev;
    cmn.io.backspace(2);
    nbufa = mbufa;
    if (kresl == 1) {
      goto statement_412;
    }
    cmn.io.backspace(4);
    read(4, fem::unformatted), eqrhs;
    cmn.io.backspace(4);
    statement_412:
    //C
    //C*** PREARE TO BACK-SUBROUTINE FROM THE CURRENT EGUATION
    //C
    ifron = npivo(nbufa);
    nikno = namev(nbufa);
    pivot = equat(ifron, nbufa);
    if (iffix(nikno) != 0) {
      vecrv(ifron) = fixed(nikno);
    }
    if (iffix(nikno) == 0) {
      equat(ifron, nbufa) = 0.0f;
    }
    //C
    //C*** BACK-SUBROUTINE IN THE CURRENT EQUATION
    //C
    FEM_DO_SAFE(jfron, 1, mfron) {
      eqrhs(nbufa) = eqrhs(nbufa) - vecrv(jfron) * equat(jfron, nbufa);
    }
    //C
    //C*** PUT THE FINAL VALUES WHERE THEY BELONG
    //C
    if (iffix(nikno) == 0) {
      vecrv(ifron) = eqrhs(nbufa) / pivot;
    }
    if (iffix(nikno) != 0) {
      fixed(nikno) = -eqrhs(nbufa);
    }
    nbufa = nbufa - 1;
    asdis(nikno) = vecrv(ifron);
  }
  //C
  //C*** ADD DISPLACEMENTS TO PREVIOUS TOTAL VALUES
  //C
  FEM_DO_SAFE(itotv, 1, ntotv) {
    tdisp(itotv) += asdis(itotv);
  }
  //C
  //C*** STORE REACTIONS FOR PRINTING LATER
  //C
  kboun = 1;
  FEM_DO_SAFE(ipoin, 1, npoin) {
    nloca = (ipoin - 1) * ndofn;
    FEM_DO_SAFE(idofn, 1, ndofn) {
      ngush = nloca + idofn;
      if (iffix(ngush) > 0) {
        goto statement_360;
      }
    }
    goto statement_370;
    statement_360:
    FEM_DO_SAFE(idofn, 1, ndofn) {
      ngash = nloca + idofn;
      treac(kboun, idofn) += fixed(ngash);
    }
    kboun++;
    statement_370:;
  }
  //C
  //C*** ADD REACTIONS INTO THE TOTAL LOAD ARRAY
  //C
  FEM_DO_SAFE(ipoin, 1, npoin) {
    FEM_DO_SAFE(ielem, 1, nelem) {
      FEM_DO_SAFE(inode, 1, nnode) {
        nloca = fem::iabs(lnods(ielem, inode));
        if (ipoin == nloca) {
          goto statement_720;
        }
      }
    }
    statement_720:
    FEM_DO_SAFE(idofn, 1, ndofn) {
      ngash = (inode - 1) * ndofn + idofn;
      mgash = (ipoin - 1) * ndofn + idofn;
      tload(ielem, ngash) += fixed(mgash);
    }
  }
}

//C$DEBUG
//C$LARGE
//C(Page 212)
void
increm(
  common& cmn,
  arr_ref<float, 2> eload,
  arr_ref<float> fixed,
  int const& iincs,
  int const& melem,
  int const& mevab,
  int& miter,
  int const& mtotv,
  int const& mvfix,
  int const& ndofn,
  int const& nelem,
  int const& nevab,
  arr_ref<int> noutp,
  arr_cref<int> nofix,
  int const& ntotv,
  int const& nvfix,
  arr_cref<float, 2> presc,
  arr_cref<float, 2> rload,
  float& tfact,
  arr_ref<float, 2> tload,
  float& toler)
{
  eload(dimension(melem, mevab));
  fixed(dimension(mtotv));
  noutp(dimension(2));
  nofix(dimension(mvfix));
  presc(dimension(mvfix, ndofn));
  rload(dimension(melem, mevab));
  tload(dimension(melem, mevab));
  common_read read(cmn);
  common_write write(cmn);
  static const char* format_900 = "('0',5x,'INCREMENT NUMBER ',i5)";
  static const char* format_960 =
    "('0',5x,'LOAD FACTOR =',f10.5,5x,' CONVERGENCE TOLERANCE=,',f10.5,/,/,5x,"
    "'MAX. NO. OF ITERATIONS =',i5,/,/,' INITIAL OUTPUT PARAMETER =',i5,5x,"
    "'FINAL OUTPUT PARAMETER =',i5)";
  //C********************************************************************
  //C
  //C*** THIS SUBROUTINE INCREMENTS THE APPLIED LOADING
  //C
  //C********************************************************************
  write(6, format_900), iincs;
  //write(6, format_900), iincs;
  //C!!
  //write(6, star), " FACTO, TOLER,MITER,NOUTP(1),NOUTP(2)";
  float facto = fem::float0;
  read(5, "(2f10.5,3i5)"), facto, toler, miter, noutp(1), noutp(2);
  //C        write(5,950) FACTO,TOLER,MITER,NOUTP(1),NOUTP(2)
  tfact += facto;
  write(6, format_960), tfact, toler, miter, noutp(1), noutp(2);
  //write(6, format_960), tfact, toler, miter, noutp(1), noutp(2);
  int ielem = fem::int0;
  int ievab = fem::int0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      eload(ielem, ievab) += rload(ielem, ievab) * facto;
      tload(ielem, ievab) += rload(ielem, ievab) * facto;
    }
  }
  //C
  //C*** INTERPRET FIXITY DATA IN VECTOR FORM
  //C
  int itotv = fem::int0;
  FEM_DO_SAFE(itotv, 1, ntotv) {
    fixed(itotv) = 0.0f;
  }
  int ivfix = fem::int0;
  int nloca = fem::int0;
  int idofn = fem::int0;
  int ngash = fem::int0;
  FEM_DO_SAFE(ivfix, 1, nvfix) {
    nloca = (nofix(ivfix) - 1) * ndofn;
    FEM_DO_SAFE(idofn, 1, ndofn) {
      ngash = nloca + idofn;
      fixed(ngash) = presc(ivfix, idofn) * facto;
    }
  }
}

//C$DEBUG
//C$LARGE
//C(Page 208)
void
input(
  common& cmn,
  arr_ref<float, 2> coord,
  arr_ref<int> iffix,
  arr_ref<int, 2> lnods,
  arr_ref<int> matno,
  int const& melem,
  int const& /* mevab */,
  int const& mfron,
  int const& mmats,
  int const& mpoin,
  int const& mtotv,
  int const& mvfix,
  int& nalgo,
  int& ncrit,
  arr_ref<int> ndfro,
  int const& ndofn,
  int& nelem,
  int& nevab,
  int& ngaus,
  int& ngau2,
  int& nincs,
  int& nmats,
  int& nnode,
  arr_ref<int> nofix,
  int& npoin,
  int const& nprop,
  int& nstre,
  int& nstr1,
  int& ntotg,
  int& ntotv,
  int& ntype,
  int& nvfix,
  arr_ref<float> posgp,
  arr_ref<float, 2> presc,
  arr_ref<float, 2> props,
  arr_ref<float> weigp)
{
  coord(dimension(mpoin, 2));
  iffix(dimension(mtotv));
  lnods(dimension(melem, 9));
  matno(dimension(melem));
  ndfro(dimension(melem));
  nofix(dimension(mvfix));
  posgp(dimension(4));
  presc(dimension(mvfix, ndofn));
  props(dimension(mmats, nprop));
  weigp(dimension(4));
  common_read read(cmn);
  common_write write(cmn);
  arr_1d<1, int> title(fem::fill0);
  int ielem = fem::int0;
  int numel = fem::int0;
  int inode = fem::int0;
  int ipoin = fem::int0;
  int idime = fem::int0;
  int ivfix = fem::int0;
  int ifpre = fem::int0;
  int idofn = fem::int0;
  int nloca = fem::int0;
  int ifdof = fem::int0;
  int ngash = fem::int0;
  int imats = fem::int0;
  int numat = fem::int0;
  int iprop = fem::int0;
  static const char* format_900 = "(11i5)";
  static const char* format_901 =
    "(/,/,' NPOIN =',i4,4x,' NELEM =',i4,4x,' NVFIX =',i4,4x,' NTYPE =',i4,4x,"
    "' NNODE =',i4,/,/,' NMATS =',i4,4x,' NGAUS =',i4,4x,' NEVAB =',i4,4x,"
    "' NALGO =',i4,/,/,' NCRIT =',i4,4x,' NINCS =',i4,4x,' NSTRE =',i4)";
  static const char* format_902 =
    "(/,/,' ELEMENT',3x,'PROPERTY',6x,'NODE NUMBERS')";
  static const char* format_903 = "(1x,i5,i9,6x,8i5)";
  static const char* format_904 = "(/,/,' NODE',10x,'X',10x,'Y')";
  static const char* format_906 = "(1x,i5,3f10.3)";
  static const char* format_907 = "(/,/,' NODE',6x,'CODE',6x,'FIXED VALUES')";
  static const char* format_908 = "(1x,i4,5x,i5,5x,5f10.6)";
  static const char* format_910 = "(/,/,' NUMBER',6x,'ELEMENT PROPERTIES')";
  static const char* format_911 = "(1x,i4,/,3x,4e14.6,/,3x,3e14.6)";
  static const char* format_920 = "(i1)";
  //C*********************************************************************
  //C
  //C**** THIS SUBROUTINE ACCEPTS MOST OF THE INPUT DATA
  //C
  //C*********************************************************************
  cmn.io.rewind(1);
  cmn.io.rewind(2);
  cmn.io.rewind(3);
  cmn.io.rewind(4);
  cmn.io.rewind(8);
  //C!!
  //write(6, star), "TITLE";
  read(5, format_920), title;
  //C        write(5,920) TITLE
  write(6, format_920), title;
  //write(6, format_920), title;
  //C
  //C*** READ THE FIRST DATA CARD, AND ECHO IT IMMEDIATELY.
  //C
  //C!!
  //write(6, star),
  //  " NPOIN,NELEM,NVFIX,NTYPE,NNODE,NMATS,NGAUS,          NALGO,NCRIT,NINCS,NS"
  //  "TRE";
  read(5, format_900), npoin, nelem, nvfix, ntype, nnode, nmats,
    ngaus, nalgo, ncrit, nincs, nstre;
  //C        write(5,900) NPOIN,NELEM,NVFIX,NTYPE,NNODE,NMATS,NGAUS,
  //C     .          NALGO,NCRIT,NINCS,NSTRE
  nevab = ndofn * nnode;
  nstr1 = nstre + 1;
  if (ntype == 3) {
    nstr1 = nstre;
  }
  ntotv = npoin * ndofn;
  ngau2 = ngaus * ngaus;
  ntotg = nelem * ngau2;
  write(6, format_901), npoin, nelem, nvfix, ntype, nnode, nmats,
    ngaus, nevab, nalgo, ncrit, nincs, nstre;
  //write(6, format_901), npoin, nelem, nvfix, ntype, nnode, nmats,
  //  ngaus, nevab, nalgo, ncrit, nincs, nstre;
  check1(cmn, ndofn, nelem, ngaus, nmats, nnode, npoin, nstre, ntype,
    nvfix, ncrit, nalgo, nincs);
  //C
  //C*** READ THE ELEMENT NODAL CONNECTIONS, AND THE PROPERTY NUMBERS.
  //C
  write(6, format_902);
  //write(6, format_902);
  FEM_DO_SAFE(ielem, 1, nelem) {
    //write(6, star), "NUMEL,MATNO(NUMEL),(LNODS(NUMEL,INODE),INODE=1,NNODE)";
    {
	  int no=0;
      read_loop rloop(cmn, 5, format_903);
      rloop, numel, no;
	  matno(numel)=no;
      FEM_DO_SAFE(inode, 1, nnode) {
        rloop, lnods(numel, inode);
      }
    }
    //C      write(5,903)NUMEL,MATNO(NUMEL),(LNODS(NUMEL,INODE),INODE=1,NNODE)
    {
      write_loop wloop(cmn, 6, format_903);
      wloop, numel, matno(numel);
      FEM_DO_SAFE(inode, 1, nnode) {
        wloop, lnods(numel, inode);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_903);
    //  wloop, numel, matno(numel);
    //  FEM_DO_SAFE(inode, 1, nnode) {
    //    wloop, lnods(numel, inode);
    //  }
    //}
  }
  //C
  //C*** ZERO ALL THE NODAL COORDINATES, PRIOR TO READING SOME OF THEM.
  //C
  FEM_DO_SAFE(ipoin, 1, npoin) {
    FEM_DO_SAFE(idime, 1, 2) {
      coord(ipoin, idime) = 0.0f;
    }
  }
  //C
  //C*** READ SOME NODAL COORDINATES, FINISHING WITH THE NODE OF ALL.
  //C
  write(6, format_904);
  //write(6, format_904);
  //C!!
  statement_6:
  //write(6, star), " IPOIN,(COORD(IPOIN,IDIME),IDIME=1,2)";
  {
    read_loop rloop(cmn, 5, "(i5,6f10.5)");
    rloop, ipoin;
    FEM_DO_SAFE(idime, 1, 2) {
      rloop, coord(ipoin, idime);
    }
  }
  //C      write(5,905) IPOIN,(COORD(IPOIN,IDIME),IDIME=1,2)
  if (ipoin != npoin) {
    goto statement_6;
  }
  //C
  //C*** INTERPOLATE COORDINATES OF MID-SIDE NODES
  //C
  nodexy(coord, lnods, melem, mpoin, nelem, nnode);
  FEM_DO_SAFE(ipoin, 1, npoin) {
    {
      write_loop wloop(cmn, 6, format_906);
      wloop, ipoin;
      FEM_DO_SAFE(idime, 1, 2) {
        wloop, coord(ipoin, idime);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_906);
    //  wloop, ipoin;
    //  FEM_DO_SAFE(idime, 1, 2) {
    //    wloop, coord(ipoin, idime);
    //  }
    //}
  }
  //C
  //C*** READ THE FIXED VALUES
  //C
  write(6, format_907);
  //write(6, format_907);
  FEM_DO_SAFE(ivfix, 1, nvfix) {
    //write(6, star), "NOFIX(IVFIX),IFPRE,(PRESC(IVFIX,IDOFN),IDOFN=1,NDOFN)";
    {
      read_loop rloop(cmn, 5, format_908);
      rloop, nofix(ivfix), ifpre;
      FEM_DO_SAFE(idofn, 1, ndofn) {
        rloop, presc(ivfix, idofn);
      }
    }
    //C      write(5,908)NOFIX(IVFIX),IFPRE,(PRESC(IVFIX,IDOFN),IDOFN=1,NDOFN)
    {
      write_loop wloop(cmn, 6, format_908);
      wloop, nofix(ivfix), ifpre;
      FEM_DO_SAFE(idofn, 1, ndofn) {
        wloop, presc(ivfix, idofn);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_908);
    //  wloop, nofix(ivfix), ifpre;
    //  FEM_DO_SAFE(idofn, 1, ndofn) {
    //    wloop, presc(ivfix, idofn);
    //  }
    //}
    nloca = (nofix(ivfix) - 1) * ndofn;
    ifdof = fem::pow(10, (ndofn - 1));
    FEM_DO_SAFE(idofn, 1, ndofn) {
      ngash = nloca + idofn;
      if (ifpre < ifdof) {
        goto statement_8;
      }
      iffix(ngash) = 1;
      ifpre = ifpre - ifdof;
      statement_8:
      ifdof = ifdof / 10;
    }
  }
  //C
  //C*** READ THE AVAILABLE SELECTION OF ELEMENT PROPERTIES.
  //C
  write(6, format_910);
  //write(6, format_910);
  FEM_DO_SAFE(imats, 1, nmats) {
    //C!!
    //write(6, star), " NUMAT";
    read(5, format_900), numat;
    //C        write(5,900) NUMAT
    //C!!
    //write(6, star), "(PROPS(NUMAT,IPROP),IPROP=1,NPROP)";
    {
      read_loop rloop(cmn, 5, "(7f10.2)");
      FEM_DO_SAFE(iprop, 1, nprop) {
        rloop, props(numat, iprop);
      }
    }
    //C        write(5,930)(PROPS(NUMAT,IPROP),IPROP=1,NPROP)
    {
      write_loop wloop(cmn, 6, format_911);
      wloop, numat;
      FEM_DO_SAFE(iprop, 1, nprop) {
        wloop, props(numat, iprop);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_911);
    //  wloop, numat;
    //  FEM_DO_SAFE(iprop, 1, nprop) {
    //    wloop, props(numat, iprop);
    //  }
    //}
  }
  //C
  //C*** SET UP GAUSSIAN INTEGRATION CONSTANTS
  //C
  gaussq(ngaus, posgp, weigp);
  check2(cmn, coord, iffix, lnods, matno, melem, mfron, mpoin, mtotv,
    mvfix, ndfro, ndofn, nelem, nmats, nnode, nofix, npoin, nvfix);
}

//C$DEBUG
//C$LARGE
//C(Page 240)
void
invar(
  arr_ref<float> devia,
  int const& lprop,
  int const& mmats,
  int const& ncrit,
  arr_cref<float, 2> props,
  float& sint3,
  float& steff,
  arr_cref<float> stemp,
  float& theta,
  float& varj2,
  float& yield)
{
  devia(dimension(4));
  props(dimension(mmats, 7));
  stemp(dimension(4));
  float root3 = fem::float0;
  float smean = fem::float0;
  float varj3 = fem::float0;
  float phira = fem::float0;
  float snphi = fem::float0;
  //C********************************************************************
  //C
  //C***THIS SUBROUTINE EVALUATED THE STRESS INVARIANTS AND THE CURRENT
  //C   VALUE OF THE YIELD FUNCTION
  //C
  //C********************************************************************
  root3 = 1.73205080757f;
  smean = (stemp(1) + stemp(2) + stemp(4)) / 3.0f;
  devia(1) = stemp(1) - smean;
  devia(2) = stemp(2) - smean;
  devia(3) = stemp(3);
  devia(4) = stemp(4) - smean;
  varj2 = devia(3) * devia(3) + 0.5f * (devia(1) * devia(1) + devia(
    2) * devia(2) + devia(4) * devia(4));
  varj3 = devia(4) * (devia(4) * devia(4) - varj2);
  steff = fem::sqrt(varj2);
  if (steff == 0.0f) {
    goto statement_10;
  }
  sint3 = -3.0f * root3 * varj3 / (2.0f * varj2 * steff);
  if (sint3 > 1.0f) {
    sint3 = 1.0f;
  }
  goto statement_20;
  statement_10:
  sint3 = 0.0f;
  statement_20:
  if (sint3 <  - 1.0f) {
    sint3 = -1.0f;
  }
  if (sint3 > 1.0f) {
    sint3 = 1.0f;
  }
  theta = fem::asin(sint3) / 3.0f;
  switch (ncrit) {
    case 1: goto statement_1;
    case 2: goto statement_2;
    case 3: goto statement_3;
    case 4: goto statement_4;
    default: break;
  }
  //C*** TRESCA
  statement_1:
  yield = 2.0f * fem::cos(theta) * steff;
  return;
  //C*** VON MISES
  statement_2:
  yield = root3 * steff;
  return;
  //C*** MOHR-COULOMB
  statement_3:
  phira = props(lprop, 7) * 0.017453292f;
  snphi = fem::sin(phira);
  yield = smean * snphi + steff * (fem::cos(theta) - fem::sin(
    theta) * snphi / root3);
  return;
  //C*** DRUCKER-PRAGER
  statement_4:
  phira = props(lprop, 7) * 0.017453292f;
  snphi = fem::sin(phira);
  //C(see pp220 (eqn.7.18)
  yield = 6.0f * smean * snphi / (root3 * (3.0f - snphi)) + steff;
}

//C$DEBUG
//C$LARGE
void
loadps(
  common& cmn,
  arr_cref<float, 2> coord,
  arr_cref<int, 2> lnods,
  arr_cref<int> matno,
  int const& melem,
  int const& mmats,
  int const& mpoin,
  int const& nelem,
  int const& nevab,
  int const& ngaus,
  int const& nnode,
  int const& npoin,
  int const& /* nstre */,
  int const& ntype,
  arr_cref<float> posgp,
  arr_cref<float, 2> props,
  arr_ref<float, 2> rload,
  arr_cref<float> weigp,
  int const& ndofn)
{
  coord(dimension(mpoin, 2));
  lnods(dimension(melem, 9));
  matno(dimension(melem));
  posgp(dimension(4));
  props(dimension(mmats, 7));
  rload(dimension(melem, 18));
  weigp(dimension(4));
  common_read read(cmn);
  common_write write(cmn);
  float twopi = fem::float0;
  int ielem = fem::int0;
  int ievab = fem::int0;
  arr_1d<1, int> title(fem::fill0);
  int iplod = fem::int0;
  int igrav = fem::int0;
  int iedge = fem::int0;
  int lodpt = fem::int0;
  arr_1d<2, float> point(fem::fill0);
  int idofn = fem::int0;
  int inode = fem::int0;
  int nloca = fem::int0;
  int ngash = fem::int0;
  float theta = fem::float0;
  float gravy = fem::float0;
  int lprop = fem::int0;
  float thick = fem::float0;
  float dense = fem::float0;
  float gxcom = fem::float0;
  float gycom = fem::float0;
  int lnode = fem::int0;
  int idime = fem::int0;
  arr_2d<2, 9, float> elcod(fem::fill0);
  int kgasp = fem::int0;
  int igaus = fem::int0;
  int jgaus = fem::int0;
  float exisp = fem::float0;
  float etasp = fem::float0;
  arr_2d<2, 9, float> deriv(fem::fill0);
  arr_1d<9, float> shape(fem::fill0);
  arr_2d<2, 9, float> cartd(fem::fill0);
  float djacb = fem::float0;
  arr_2d<2, 9, float> gpcod(fem::fill0);
  float dvolu = fem::float0;
  int mgash = fem::int0;
  int nedge = fem::int0;
  int nodeg = fem::int0;
  int ncode = fem::int0;
  int neass = fem::int0;
  arr_1d<4, int> noprs(fem::fill0);
  int iodeg = fem::int0;
  arr_2d<4, 2, float> press(fem::fill0);
  arr_1d<2, float> pgash(fem::fill0);
  arr_1d<2, float> dgash(fem::fill0);
  float pxcom = fem::float0;
  float pycom = fem::float0;
  float radus = fem::float0;
  int jnode = fem::int0;
  int kount = fem::int0;
  int knode = fem::int0;
  static const char* format_905 = "(1x,i3,2x,8f9.2,/(5x,8f9.2))";
  static const char* format_907 =
    "('0',5x,' TOTAL NODAL FORCES FOR EACH ELEMENT')";
  static const char* format_911 =
    "('0','GRAVITY ANGLE =',f10.3,' GRAVITY CONSTANT =',f10.3)";
  static const char* format_912 = "('0',5x,'NO. OF LOADED EDGES =',i5)";
  static const char* format_913 = "(i10,5x,3i5)";
  static const char* format_914 = "(6f10.3)";
  static const char* format_915 =
    "('0',5x,'LIST OF LOADED EDGES AND APPLIED LOADS')";
  static const char* format_919 =
    "(1x,'IPLOD=',i5,5x,'IGRAV=',i5,5x,'IEDGE=',i5)";
  static const char* format_931 = "(i5,2f10.3)";
  //C********************************************************************
  //C
  //C****THIS SUBROUTINE EVALUATES THE CONSISTENT NODAL FORCES FOR EACH
  //C    ELEMENT
  //C
  //C********************************************************************
  twopi = 6.283185308f;
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      rload(ielem, ievab) = 0.0f;
    }
  }
  //C!!
  //write(6, star), " TITLE";
  read(5, "(2x,i1)"), title;
  //C        write(5,901) TITLE
  write(6, "(2x,i1)"), title;
  //C
  //C*** READ DATA CONTROLLING LOADING TYPES TO BE INPUTTED
  //C
  //C!!
  //write(6, star), " IPLOD,IGRAV,IEDGE";
  read(5, "(3i5)"), iplod, igrav, iedge;
  //C        write(5,918) IPLOD,IGRAV,IEDGE
  write(6, format_919), iplod, igrav, iedge;
  //write(6, format_919), iplod, igrav, iedge;
  //C
  //C*** READ NODAL POINT LOADS
  //C
  if (iplod == 0) {
    goto statement_500;
  }
  //C!!
  statement_20:
  //write(6, star), "LODPT,(POINT(IDOFN),IDOFN=1,2)";
  {
    read_loop rloop(cmn, 5, format_931);
    rloop, lodpt;
    FEM_DO_SAFE(idofn, 1, 2) {
      rloop, point(idofn);
    }
  }
  //C      write(5,931)LODPT,(POINT(IDOFN),IDOFN=1,2)
  {
    write_loop wloop(cmn, 6, format_931);
    wloop, lodpt;
    FEM_DO_SAFE(idofn, 1, 2) {
      wloop, point(idofn);
    }
  }
  //{
  //  write_loop wloop(cmn, 6, format_931);
  //  wloop, lodpt;
  //  FEM_DO_SAFE(idofn, 1, 2) {
  //    wloop, point(idofn);
  //  }
  //}
  //C
  //C*** ASSOCIATE THE NODAL POINT LOADS WITH AN ELEMENT
  //C
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(inode, 1, nnode) {
      nloca = fem::iabs(lnods(ielem, inode));
      if (lodpt == nloca) {
        goto statement_40;
      }
    }
  }
  statement_40:
  FEM_DO_SAFE(idofn, 1, 2) {
    ngash = (inode - 1) * 2 + idofn;
    rload(ielem, ngash) = point(idofn);
  }
  if (lodpt < npoin) {
    goto statement_20;
  }
  statement_500:
  if (igrav == 0) {
    goto statement_600;
  }
  //C
  //C*** GRAVITY LOADING SECSION
  //C
  //C*** READ GRAVITY ANGLE AND GRAVITATIONAL CONSTANT
  //C
  //C!!
  //write(6, star), "THETA,GRAVY";
  //read(6, star), theta, gravy;
  write(5, "(2f10.3)"), theta, gravy;
  write(6, format_911), theta, gravy;
  //write(6, format_911), theta, gravy;
  theta = theta / 57.295779514f;
  //C
  //C*** LOOP OVER EACH ELEMENT
  //C
  FEM_DO_SAFE(ielem, 1, nelem) {
    //C
    //C*** SET UP PRELIMINARY CONSTANTS
    //C
    lprop = matno(ielem);
    thick = props(lprop, 3);
    dense = props(lprop, 4);
    if (dense == 0.0f) {
      goto statement_90;
    }
    gxcom = dense * gravy * fem::sin(theta);
    gycom = -dense * gravy * fem::cos(theta);
    //C
    //C*** COMPUTE COORDINATES OF THE ELEMENT NODAL POINTS
    //C
    FEM_DO_SAFE(inode, 1, nnode) {
      lnode = fem::iabs(lnods(ielem, inode));
      FEM_DO_SAFE(idime, 1, 2) {
        elcod(idime, inode) = coord(lnode, idime);
      }
    }
    //C
    //C*** ENTER LOOPS FOR AREA NUMERICAL INTEGRATION
    //C
    kgasp = 0;
    FEM_DO_SAFE(igaus, 1, ngaus) {
      FEM_DO_SAFE(jgaus, 1, ngaus) {
        exisp = posgp(igaus);
        etasp = posgp(jgaus);
        //C
        //C*** COMPUTE THE SHAPE FUNCTIONS AT THE SAMPLING POINTS AND ELEMENTAL
        //C    VOLUME
        //C
        sfr2(deriv, etasp, exisp, nnode, shape);
        kgasp++;
        jacob2(cmn, cartd, deriv, djacb, elcod, gpcod, ielem, kgasp,
          nnode, shape);
        dvolu = djacb * weigp(igaus) * weigp(jgaus);
        if (thick != 0.0f) {
          dvolu = dvolu * thick;
        }
        if (ntype == 3) {
          dvolu = dvolu * twopi * gpcod(1, kgasp);
        }
        //C
        //C*** CALCULATE LOADS AND ASSOCIAT WITH ELEMENT NODAL POINTS
        //C
        FEM_DO_SAFE(inode, 1, nnode) {
          ngash = (inode - 1) * 2 + 1;
          mgash = (inode - 1) * 2 + 2;
          rload(ielem, ngash) += gxcom * shape(inode) * dvolu;
          rload(ielem, mgash) += gycom * shape(inode) * dvolu;
        }
      }
    }
    statement_90:;
  }
  statement_600:
  if (iedge == 0) {
    goto statement_700;
  }
  //C
  //C*** DISTRIBUTED EDGE LOADS SECTION
  //C
  //C!!
  //write(6, star), "NEDGE";
  read(5, "(i5)"), nedge;
  //C        write(5,932)NEDGE
  write(6, format_912), nedge;
  //write(6, format_912), nedge;
  write(6, format_915);
  //write(6, format_915);
  nodeg = 3;
  ncode = nnode;
  if (nnode == 4) {
    nodeg = 2;
  }
  if (nnode == 9) {
    ncode = 8;
  }
  //C
  //C*** LOOP OVER EACH LOADED EDGE
  //C
  FEM_DO_SAFE(iedge, 1, nedge) {
    //C
    //C*** READ DATA LOCATING THE LOADED EDGE AND APPLIED LOAD
    //C!!
    //write(6, star), " NEASS,(NOPRS(IODEG),IODEG=1,NODEG)";
    {
      read_loop rloop(cmn, 5, "(4i5)");
      rloop, neass;
      FEM_DO_SAFE(iodeg, 1, nodeg) {
        rloop, noprs(iodeg);
      }
    }
    //C        write(5,902) NEASS,(NOPRS(IODEG),IODEG=1,NODEG)
    {
      write_loop wloop(cmn, 6, format_913);
      wloop, neass;
      FEM_DO_SAFE(iodeg, 1, nodeg) {
        wloop, noprs(iodeg);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_913);
    //  wloop, neass;
    //  FEM_DO_SAFE(iodeg, 1, nodeg) {
    //    wloop, noprs(iodeg);
    //  }
    //}
    //C!!
    //write(6, star), "((PRESS(IODEG,IDOFN),IDOFN=1,2),IODEG=1,NODEG)";
    {
      read_loop rloop(cmn, 5, format_914);
      FEM_DO_SAFE(iodeg, 1, nodeg) {
        FEM_DO_SAFE(idofn, 1, 2) {
          rloop, press(iodeg, idofn);
        }
      }
    }
    //C        write(5,914)((PRESS(IODEG,IDOFN),IDOFN=1,2),IODEG=1,NODEG)
    {
      write_loop wloop(cmn, 6, format_914);
      FEM_DO_SAFE(iodeg, 1, nodeg) {
        FEM_DO_SAFE(idofn, 1, 2) {
          wloop, press(iodeg, idofn);
        }
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_914);
    //  FEM_DO_SAFE(iodeg, 1, nodeg) {
    //    FEM_DO_SAFE(idofn, 1, 2) {
    //      wloop, press(iodeg, idofn);
    //    }
    //  }
    //}
    etasp = -1.0f;
    //C
    //C*** CALCULITE THE COORDINATES OF THE NODES OF THE ELEMENT EDGE
    //C
    FEM_DO_SAFE(iodeg, 1, nodeg) {
      lnode = noprs(iodeg);
      FEM_DO_SAFE(idime, 1, 2) {
        elcod(idime, iodeg) = coord(lnode, idime);
      }
    }
    //C
    //C*** ENTER LOOP FOR LINEAR NUMERICAL INTEGRATION
    //C
    FEM_DO_SAFE(igaus, 1, ngaus) {
      exisp = posgp(igaus);
      //C
      //C*** EVALUATE THE SHAPE FUNCTIOS AT THE SAMPLING POINTS
      //C
      sfr2(deriv, etasp, exisp, nnode, shape);
      //C
      //C*** CALCULATE COMPNENTS OF THE EQUIVALENT NODAL LOADS
      //C
      FEM_DO_SAFE(idofn, 1, 2) {
        pgash(idofn) = 0.0f;
        dgash(idofn) = 0.0f;
        FEM_DO_SAFE(iodeg, 1, nodeg) {
          pgash(idofn) += press(iodeg, idofn) * shape(iodeg);
          dgash(idofn) += elcod(idofn, iodeg) * deriv(1, iodeg);
        }
      }
      dvolu = weigp(igaus);
      pxcom = dgash(1) * pgash(2) - dgash(2) * pgash(1);
      pycom = dgash(1) * pgash(1) + dgash(2) * pgash(2);
      if (ntype != 3) {
        goto statement_115;
      }
      radus = 0.0f;
      FEM_DO_SAFE(iodeg, 1, nodeg) {
        radus += shape(iodeg) * elcod(1, iodeg);
      }
      dvolu = dvolu * twopi * radus;
      statement_115:
      //C
      //C*** ASSOCIATE THE EQUIVALENT NODAL EDGE LOADS WITH AN ELEMENT
      //C
      FEM_DO_SAFE(inode, 1, nnode) {
        nloca = fem::iabs(lnods(neass, inode));
        if (nloca == noprs(1)) {
          goto statement_130;
        }
      }
      statement_130:
      jnode = inode + nodeg - 1;
      kount = 0;
      FEM_DO_SAFE(knode, inode, jnode) {
        kount++;
        ngash = (knode - 1) * ndofn + 1;
        mgash = (knode - 1) * ndofn + 2;
        if (knode > ncode) {
          ngash = 1;
        }
        if (knode > ncode) {
          mgash = 2;
        }
        rload(neass, ngash) += shape(kount) * pxcom * dvolu;
        rload(neass, mgash) += shape(kount) * pycom * dvolu;
      }
    }
  }
  statement_700:
  write(6, format_907);
  //write(6, format_907);
  FEM_DO_SAFE(ielem, 1, nelem) {
    {
      write_loop wloop(cmn, 6, format_905);
      wloop, ielem;
      FEM_DO_SAFE(ievab, 1, nevab) {
        wloop, rload(ielem, ievab);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_905);
    //  wloop, ielem;
    //  FEM_DO_SAFE(ievab, 1, nevab) {
    //    wloop, rload(ielem, ievab);
    //  }
    //}
  }
}

//C$DEBUG
//C$LARGE
//C(Page 258)
void
output(
  common& cmn,
  int const& iiter,
  int const& mtotg,
  int const& mtotv,
  int const& mvfix,
  int const& nelem,
  int const& ngaus,
  arr_cref<int> nofix,
  arr_cref<int> noutp,
  int const& npoin,
  int const& nvfix,
  arr_cref<float, 2> strsg,
  arr_cref<float> tdisp,
  arr_cref<float, 2> treac,
  arr_cref<float> epstn,
  int const& ntype,
  int const& nchek)
{
  nofix(dimension(mvfix));
  noutp(dimension(2));
  strsg(dimension(4, mtotg));
  tdisp(dimension(mtotv));
  treac(dimension(mvfix, 2));
  epstn(dimension(mtotg));
  common_write write(cmn);
  int koutp = fem::int0;
  int ipoin = fem::int0;
  int ngash = fem::int0;
  int ngish = fem::int0;
  int igash = fem::int0;
  int ivfix = fem::int0;
  int idofn = fem::int0;
  int kgaus = fem::int0;
  int ielem = fem::int0;
  int kelgs = fem::int0;
  int igaus = fem::int0;
  int jgaus = fem::int0;
  float xgash = fem::float0;
  float xgish = fem::float0;
  float xgesh = fem::float0;
  float xgosh = fem::float0;
  arr_1d<3, float> strsp(fem::fill0);
  int istr1 = fem::int0;
  int istre = fem::int0;
  static const char* format_900 = "('0',6x,'DISPLACEMENTS')";
  static const char* format_910 = "(i10,3e14.6)";
  static const char* format_920 = "('0',5x,'REACTIONS')";
  static const char* format_930 = "('0',5x,'ELEMENT NO. =',i5)";
  static const char* format_940 = "(i5,2x,6e14.6,f8.3,e14.6)";
  static const char* format_950 = "('0',6x,'NODE',6x,'X-DISP.',7x,'Y-DISP.')";
  static const char* format_970 =
    "('0',1x,'G.P.',6x,'XX-STRESS',5x,'YY-STRESS',5x,'XY-STRESS',5x,"
    "'ZZ-STRESS',6x,'MAX P.S.',6x,'MIN P.S.',3x,'ANGLE',3x,'E.P.S.')";
  static const char* format_975 =
    "('0',1x,'G.P.',6x,'RR-STRESS',5x,'ZZ-STRESS',5x,'RZ-STRESS',5x,"
    "'TT-STRESS',6x,'MAX P.S.',6x,'MIN P.S.',3x,'ANGLE',3x,'E.P.S.')";
  //C*******************************************************************
  //C
  //C**** THIS SUBROUTINE OUTPUT DISPLACEMENTS.REACTION AND STRESSES
  //C
  //C*******************************************************************
  koutp = noutp(1);
  if (iiter > 1) {
    koutp = noutp(2);
  }
  if (iiter == 1 && nchek == 0) {
    koutp = noutp(2);
  }
  //C
  //C*** OUTPUT DISPLACEMENTS
  //C
  if (koutp < 1) {
    goto statement_10;
  }
  write(6, format_900);
  //write(6, format_900);
  //if (ntype != 3) {
  //  write(6, format_950);
  //}
  if (ntype != 3) {
    write(6, format_950);
  }
  if (ntype == 3) {
    write(6, "('0',6x,'NODE',6x,'R-DISP.',7x,'Z-DISP.')");
  }
  FEM_DO_SAFE(ipoin, 1, npoin) {
    ngash = ipoin * 2;
    ngish = ngash - 2 + 1;
    {
      write_loop wloop(cmn, 6, format_910);
      wloop, ipoin;
      FEM_DO_SAFE(igash, ngish, ngash) {
        wloop, tdisp(igash);
      }
    }
    //{
    //  write_loop wloop(cmn, 6, format_910);
    //  wloop, ipoin;
    //  FEM_DO_SAFE(igash, ngish, ngash) {
    //    wloop, tdisp(igash);
    //  }
    //}
  }
  statement_10:
  //C
  //C*** OUTPUT REACTIONS
  //C
  if (koutp < 2) {
    goto statement_30;
  }
  write(6, format_920);
  //write(6, format_920);
  if (ntype != 3) {
    write(6, "('0',6x,'NODE',6x,'X-REAC.',7x,'Y-REAC.')");
  }
  if (ntype == 3) {
    write(6, "('0',6x,'NODE',6x,'R-REAC.',7x,'Z-REAC.')");
  }
  FEM_DO_SAFE(ivfix, 1, nvfix) {
    {
      write_loop wloop(cmn, 6, format_910);
      wloop, nofix(ivfix);
      FEM_DO_SAFE(idofn, 1, 2) {
        wloop, treac(ivfix, idofn);
      }
    }
  }
  statement_30:
  //C
  //C*** OUTPUT
  //C
  if (koutp < 3) {
    goto statement_50;
  }
  if (ntype != 3) {
    write(6, format_970);
  }
  //if (ntype != 3) {
  //  write(6, format_970);
  //}
  if (ntype == 3) {
    write(6, format_975);
  }
  //if (ntype == 3) {
  //  write(6, format_975);
  //}
  kgaus = 0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    kelgs = 0;
    write(6, format_930), ielem;
    //write(6, format_930), ielem;
    FEM_DO_SAFE(igaus, 1, ngaus) {
      FEM_DO_SAFE(jgaus, 1, ngaus) {
        kgaus++;
        kelgs++;
        xgash = (strsg(1, kgaus) + strsg(2, kgaus)) * 0.5f;
        xgish = (strsg(1, kgaus) - strsg(2, kgaus)) * 0.5f;
        xgesh = strsg(3, kgaus);
        xgosh = fem::sqrt(xgish * xgish + xgesh * xgesh);
        strsp(1) = xgash + xgosh;
        strsp(2) = xgash - xgosh;
        if (xgish == 0.0f) {
          xgish = 0.1e-20f;
        }
        strsp(3) = fem::atan(xgesh / xgish) * 28.647889757f;
        {
          write_loop wloop(cmn, 6, format_940);
          wloop, kelgs;
          FEM_DO_SAFE(istr1, 1, 4) {
            wloop, strsg(istr1, kgaus);
          }
          FEM_DO_SAFE(istre, 1, 3) {
            wloop, strsp(istre);
          }
          wloop, epstn(kgaus);
        }
        //{
        //  write_loop wloop(cmn, 6, format_940);
        //  wloop, kelgs;
        //  FEM_DO_SAFE(istr1, 1, 4) {
        //    wloop, strsg(istr1, kgaus);
        //  }
        //  FEM_DO_SAFE(istre, 1, 3) {
        //    wloop, strsp(istre);
        //  }
        //  wloop, epstn(kgaus);
        //}
      }
    }
  }
  statement_50:;
}

//C$DEBUG
//C$LARGE
//C(Page 2
void
linear(
  arr_cref<float, 2> cartd,
  arr_cref<float, 2> dmatx,
  arr_cref<float, 2> eldis,
  int const& lprop,
  int const& mmats,
  int const& ndofn,
  int const& nnode,
  int const& nstre,
  int const& ntype,
  arr_cref<float, 2> props,
  arr_ref<float> stran,
  arr_ref<float> stres,
  int const& kgasp,
  arr_cref<float, 2> gpcod,
  arr_cref<float> shape)
{
  cartd(dimension(2, 9));
  dmatx(dimension(4, 4));
  eldis(dimension(2, 9));
  props(dimension(mmats, 7));
  stran(dimension(4));
  stres(dimension(4));
  gpcod(dimension(2, 9));
  shape(dimension(9));
  //C*******************************************************************
  //C
  //C**** THIS SUBROUTINE EVALUATES STRESSES AND STRAINS ASSNMING LINEAR
  //C     ELASTIC BEHAVIOUR
  //C
  //C*******************************************************************
  float poiss = props(lprop, 2);
  int idofn = fem::int0;
  int jdofn = fem::int0;
  float bgash = fem::float0;
  int inode = fem::int0;
  arr_2d<2, 2, float> agash(fem::fill0);
  FEM_DO_SAFE(idofn, 1, ndofn) {
    FEM_DO_SAFE(jdofn, 1, ndofn) {
      bgash = 0.0f;
      FEM_DO_SAFE(inode, 1, nnode) {
        bgash += cartd(jdofn, inode) * eldis(idofn, inode);
      }
      agash(idofn, jdofn) = bgash;
    }
  }
  //C
  //C*** CALCULATE THE STRAINS
  //C
  stran(1) = agash(1, 1);
  stran(2) = agash(2, 2);
  stran(3) = agash(1, 2) + agash(2, 1);
  stran(4) = 0.0f;
  FEM_DO_SAFE(inode, 1, nnode) {
    stran(4) += eldis(1, inode) * shape(inode) / gpcod(1, kgasp);
  }
  //C
  //C*** AND THE CORRESPONDING STRESSES
  //C
  int istre = fem::int0;
  int jstre = fem::int0;
  FEM_DO_SAFE(istre, 1, nstre) {
    stres(istre) = 0.0f;
    FEM_DO_SAFE(jstre, 1, nstre) {
      stres(istre) += dmatx(istre, jstre) * stran(jstre);
    }
  }
  if (ntype == 1) {
    stres(4) = 0.0f;
  }
  if (ntype == 2) {
    stres(4) = poiss * (stres(1) + stres(2));
  }
}

//C$DEBUG
//C$LARGE
//C(Page 241)
void
yieldf(
  arr_ref<float> avect,
  arr_cref<float> devia,
  int const& lprop,
  int const& mmats,
  int const& ncrit,
  int const& nstr1,
  arr_cref<float, 2> props,
  float const& /* sint3 */,
  float const& steff,
  float const& theta,
  float const& varj2)
{
  avect(dimension(4));
  devia(dimension(4));
  props(dimension(mmats, 7));
  float frict = fem::float0;
  float tanth = fem::float0;
  float tant3 = fem::float0;
  float sinth = fem::float0;
  float costh = fem::float0;
  float cost3 = fem::float0;
  float root3 = fem::float0;
  arr_1d<4, float> veca1(fem::fill0);
  int istr1 = fem::int0;
  arr_1d<4, float> veca2(fem::fill0);
  arr_1d<4, float> veca3(fem::fill0);
  float cons1 = fem::float0;
  float abthe = fem::float0;
  float cons2 = fem::float0;
  float cons3 = fem::float0;
  float plumi = fem::float0;
  float snphi = fem::float0;
  //C******************************************************************
  //C
  //C**** THIS SUBROUTINE EVALUATES THE FLOW VECTOR
  //C
  //C******************************************************************
  if (steff == 0.0f) {
    return;
  }
  frict = props(lprop, 7);
  tanth = fem::tan(theta);
  tant3 = fem::tan(3.0f * theta);
  sinth = fem::sin(theta);
  costh = fem::cos(theta);
  cost3 = fem::cos(3.0f * theta);
  root3 = 1.73205080757f;
  //C
  //C*** CALCULATE VECTOR A1
  //C
  veca1(1) = 1.0f;
  veca1(2) = 1.0f;
  veca1(3) = 0.0f;
  veca1(4) = 1.0f;
  //C
  //C*** CALCULATE VECTOR A2
  //C
  FEM_DO_SAFE(istr1, 1, nstr1) {
    veca2(istr1) = devia(istr1) / (2.0f * steff);
  }
  veca2(3) = devia(3) / steff;
  //C
  //C*** CALCULATE VECTOR A3
  //C
  veca3(1) = devia(2) * devia(4) + varj2 / 3.0f;
  veca3(2) = devia(1) * devia(4) + varj2 / 3.0f;
  veca3(3) = -2.0f * devia(3) * devia(4);
  veca3(4) = devia(1) * devia(2) - devia(3) * devia(3) + varj2 / 3.0f;
  switch (ncrit) {
    case 1: goto statement_1;
    case 2: goto statement_2;
    case 3: goto statement_3;
    case 4: goto statement_4;
    default: break;
  }
  //C
  //C*** TRESCA
  //C
  statement_1:
  cons1 = 0.0f;
  abthe = fem::abs(theta * 57.29577951308f);
  if (abthe < 29.0f) {
    goto statement_20;
  }
  cons2 = root3;
  cons3 = 0.0f;
  goto statement_40;
  statement_20:
  cons2 = 2.0f * (costh + sinth * tant3);
  cons3 = root3 * sinth / (varj2 * cost3);
  goto statement_40;
  //C
  //C*** VON MISES
  //C
  statement_2:
  cons1 = 0.0f;
  cons2 = root3;
  cons3 = 0.0f;
  goto statement_40;
  //C
  //C*** MOHR-COULOMB
  //C
  statement_3:
  cons1 = fem::sin(frict * 0.017453292f) / 3.0f;
  abthe = fem::abs(theta * 57.29577951308f);
  if (abthe < 29.0f) {
    goto statement_30;
  }
  cons3 = 0.0f;
  plumi = 1.0f;
  if (theta > 0.0f) {
    plumi = -1.0f;
  }
  cons2 = 0.5f * (root3 + plumi * cons1 * root3);
  goto statement_40;
  statement_30:
  cons2 = costh * ((1.0f + tanth * tant3) + cons1 * (tant3 - tanth) * root3);
  cons3 = (root3 * sinth + 3.0f * cons1 * costh) / (2.0f * varj2 * cost3);
  goto statement_40;
  //C
  //C*** DRUCKER-PRAGER
  //C
  statement_4:
  snphi = fem::sin(frict * 0.017453292f);
  cons1 = 2.0f * snphi / (root3 * (3.0f - snphi));
  cons2 = 1.0f;
  cons3 = 0.0f;
  statement_40:
  FEM_DO_SAFE(istr1, 1, nstr1) {
    avect(istr1) = cons1 * veca1(istr1) + cons2 * veca2(istr1) +
      cons3 * veca3(istr1);
  }
}

//C
//C$DEBUG
//C$LARGE
//C(Page 253)
void
residu(
  common& cmn,
  arr_cref<float> asdis,
  arr_cref<float, 2> coord,
  arr_ref<float> effst,
  arr_ref<float, 2> eload,
  float const& /* facto */,
  int const& /* iiter */,
  arr_cref<int, 2> lnods,
  int& lprop,
  arr_cref<int> matno,
  int const& melem,
  int const& mmats,
  int const& mpoin,
  int const& mtotg,
  int const& mtotv,
  int const& ndofn,
  int const& nelem,
  int const& nevab,
  int const& ngaus,
  int const& nnode,
  int const& nstr1,
  int const& ntype,
  arr_cref<float> posgp,
  arr_cref<float, 2> props,
  int const& nstre,
  int const& ncrit,
  arr_ref<float, 2> strsg,
  arr_cref<float> weigp,
  arr_cref<float> /* tdisp */,
  arr_ref<float> epstn)
{
  asdis(dimension(mtotv));
  coord(dimension(mpoin, 2));
  effst(dimension(mtotg));
  eload(dimension(melem, 18));
  lnods(dimension(melem, 9));
  matno(dimension(melem));
  posgp(dimension(4));
  props(dimension(mmats, 7));
  strsg(dimension(4, mtotg));
  weigp(dimension(4));
  epstn(dimension(mtotg));
  float root3 = fem::float0;
  float twopi = fem::float0;
  int ielem = fem::int0;
  int ievab = fem::int0;
  int kgaus = fem::int0;
  float uniax = fem::float0;
  float hards = fem::float0;
  float frict = fem::float0;
  int inode = fem::int0;
  int lnode = fem::int0;
  int nposn = fem::int0;
  int idofn = fem::int0;
  arr_2d<2, 9, float> elcod(fem::fill0);
  arr_2d<2, 9, float> eldis(fem::fill0);
  arr_2d<4, 4, float> dmatx(fem::fill0);
  float thick = fem::float0;
  int kgasp = fem::int0;
  int igaus = fem::int0;
  int jgaus = fem::int0;
  float exisp = fem::float0;
  float etasp = fem::float0;
  arr_2d<2, 9, float> deriv(fem::fill0);
  arr_1d<9, float> shape(fem::fill0);
  arr_2d<2, 9, float> cartd(fem::fill0);
  float djacb = fem::float0;
  arr_2d<2, 9, float> gpcod(fem::fill0);
  float dvolu = fem::float0;
  arr_2d<4, 18, float> bmatx(fem::fill0);
  arr_1d<4, float> stran(fem::fill0);
  arr_1d<4, float> stres(fem::fill0);
  float preys = fem::float0;
  int istr1 = fem::int0;
  arr_1d<4, float> desig(fem::fill0);
  arr_1d<4, float> sigma(fem::fill0);
  arr_1d<4, float> devia(fem::fill0);
  float sint3 = fem::float0;
  float steff = fem::float0;
  float theta = fem::float0;
  float varj2 = fem::float0;
  float yield = fem::float0;
  float espre = fem::float0;
  float escur = fem::float0;
  float rfact = fem::float0;
  int mstep = fem::int0;
  float astep = fem::float0;
  float reduc = fem::float0;
  arr_1d<4, float> sgtot(fem::fill0);
  int istep = fem::int0;
  arr_1d<4, float> avect(fem::fill0);
  float abeta = fem::float0;
  arr_1d<4, float> dvect(fem::fill0);
  float agash = fem::float0;
  float dlamd = fem::float0;
  float bgash = fem::float0;
  float curys = fem::float0;
  float bring = fem::float0;
  int mgash = fem::int0;
  int istre = fem::int0;
  //C********************************************************************
  //C
  //C**** THIS SUBROUTINE REDUCES THE STRESSES TO THE YIELD SURFACE AND
  //C      EVALUATES THE EQUIVALENT NODAL FORCES
  //C
  //C**********************************************************************
  root3 = 1.73205080757f;
  twopi = 6.283185308f;
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      eload(ielem, ievab) = 0.0f;
    }
  }
  kgaus = 0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    lprop = matno(ielem);
    uniax = props(lprop, 5);
    hards = props(lprop, 6);
    frict = props(lprop, 7);
    if (ncrit == 3) {
      uniax = props(lprop, 5) * fem::cos(frict * 0.017453292f);
    }
    if (ncrit == 4) {
      uniax = 6.0f * props(lprop, 5) * fem::cos(frict *
        0.017453292f) / (root3 * (3.0f - fem::sin(frict *
        0.017453292f)));
    }
    //C
    //C*** COMPUTE COORDINATE AND INCREMENTAL DISPLACEMENTS OF THE
    //C    ELEMENT NODAL POINTS
    //C
    FEM_DO_SAFE(inode, 1, nnode) {
      lnode = fem::iabs(lnods(ielem, inode));
      nposn = (lnode - 1) * ndofn;
      FEM_DO_SAFE(idofn, 1, ndofn) {
        nposn++;
        elcod(idofn, inode) = coord(lnode, idofn);
        eldis(idofn, inode) = asdis(nposn);
      }
    }
    modps(dmatx, lprop, mmats, ntype, props);
    thick = props(lprop, 3);
    kgasp = 0;
    FEM_DO_SAFE(igaus, 1, ngaus) {
      FEM_DO_SAFE(jgaus, 1, ngaus) {
        exisp = posgp(igaus);
        etasp = posgp(jgaus);
        kgaus++;
        kgasp++;
        sfr2(deriv, etasp, exisp, nnode, shape);
        jacob2(cmn, cartd, deriv, djacb, elcod, gpcod, ielem, kgasp,
          nnode, shape);
        dvolu = djacb * weigp(igaus) * weigp(jgaus);
        if (ntype == 3) {
          dvolu = dvolu * twopi * gpcod(1, kgasp);
        }
        if (thick != 0.0f) {
          dvolu = dvolu * thick;
        }
        bmatps(bmatx, cartd, nnode, shape, gpcod, ntype, kgasp);
        linear(cartd, dmatx, eldis, lprop, mmats, ndofn, nnode,
          nstre, ntype, props, stran, stres, kgasp, gpcod, shape);
        preys = uniax + epstn(kgaus) * hards;
        FEM_DO_SAFE(istr1, 1, nstr1) {
          desig(istr1) = stres(istr1);
          sigma(istr1) = strsg(istr1, kgaus) + stres(istr1);
        }
        invar(devia, lprop, mmats, ncrit, props, sint3, steff, sigma,
          theta, varj2, yield);
        espre = effst(kgaus) - preys;
        if (espre >= 0.0f) {
          goto statement_50;
        }
        escur = yield - preys;
        if (escur <= 0.0f) {
          goto statement_60;
        }
        rfact = escur / (yield - effst(kgaus));
        goto statement_70;
        statement_50:
        escur = yield - effst(kgaus);
        if (escur <= 0.0f) {
          goto statement_60;
        }
        rfact = 1.0f;
        statement_70:
        mstep = escur * 8.0f / uniax + 1.0f;
        astep = mstep;
        reduc = 1.0f - rfact;
        FEM_DO_SAFE(istr1, 1, nstr1) {
          sgtot(istr1) = strsg(istr1, kgaus) + reduc * stres(istr1);
          stres(istr1) = rfact * stres(istr1) / astep;
        }
        FEM_DO_SAFE(istep, 1, mstep) {
          invar(devia, lprop, mmats, ncrit, props, sint3, steff,
            sgtot, theta, varj2, yield);
          yieldf(avect, devia, lprop, mmats, ncrit, nstr1, props,
            sint3, steff, theta, varj2);
          flowpl(avect, abeta, dvect, ntype, props, lprop, nstr1, mmats);
          agash = 0.0f;
          FEM_DO_SAFE(istr1, 1, nstr1) {
            agash += avect(istr1) * stres(istr1);
          }
          dlamd = agash * abeta;
          if (dlamd < 0.0f) {
            dlamd = 0.0f;
          }
          bgash = 0.0f;
          FEM_DO_SAFE(istr1, 1, nstr1) {
            bgash += avect(istr1) * sgtot(istr1);
            sgtot(istr1) += stres(istr1) - dlamd * dvect(istr1);
          }
          epstn(kgaus) += dlamd * bgash / yield;
        }
        invar(devia, lprop, mmats, ncrit, props, sint3, steff, sgtot,
          theta, varj2, yield);
        curys = uniax + epstn(kgaus) * hards;
        bring = 1.0f;
        if (yield > curys) {
          bring = curys / yield;
        }
        FEM_DO_SAFE(istr1, 1, nstr1) {
          strsg(istr1, kgaus) = bring * sgtot(istr1);
        }
        effst(kgaus) = bring * yield;
        //C*** ALTERNATIVE LOCATION OF STRESS REDUCTION LOOP TERMINATION CARD
        //C  90 CONTINUE
        //C***
        goto statement_190;
        statement_60:
        FEM_DO_SAFE(istr1, 1, nstr1) {
          strsg(istr1, kgaus) += desig(istr1);
        }
        effst(kgaus) = yield;
        //C
        //C*** CALCULATE THE EQUIVALENT NODAL FOCES AND ASSOCIATE WITH THE
        //C    ELEMENT NODES
        statement_190:
        mgash = 0;
        FEM_DO_SAFE(inode, 1, nnode) {
          FEM_DO_SAFE(idofn, 1, ndofn) {
            mgash++;
            FEM_DO_SAFE(istre, 1, nstre) {
              eload(ielem, mgash) += bmatx(istre, mgash) * strsg(istre,
                kgaus) * dvolu;
            }
          }
        }
      }
    }
  }
}

//C$debug
//C$large
//C(Pa
void
stiffp(
  common& cmn,
  arr_cref<float, 2> coord,
  arr_cref<float> epstn,
  int const& iincs,
  arr_cref<int, 2> lnods,
  arr_cref<int> matno,
  int const& mevab,
  int const& mmats,
  int const& mpoin,
  int const& /* mtotv */,
  int const& nelem,
  int const& nevab,
  int const& ngaus,
  int const& nnode,
  int const& nstre,
  int& nstr1,
  arr_cref<float> posgp,
  arr_cref<float, 2> props,
  arr_cref<float> weigp,
  int const& melem,
  int const& mtotg,
  arr_cref<float, 2> strsg,
  int const& ntype,
  int const& ncrit)
{
  coord(dimension(mpoin, 2));
  epstn(dimension(mtotg));
  lnods(dimension(melem, 9));
  matno(dimension(melem));
  posgp(dimension(4));
  props(dimension(mmats, 7));
  weigp(dimension(4));
  strsg(dimension(4, mtotg));
  common_write write(cmn);
  float twopi = fem::float0;
  int kgaus = fem::int0;
  int ielem = fem::int0;
  int lprop = fem::int0;
  int inode = fem::int0;
  int lnode = fem::int0;
  int iposn = fem::int0;
  int idime = fem::int0;
  arr_2d<2, 9, float> elcod(fem::fill0);
  float thick = fem::float0;
  int ievab = fem::int0;
  int jevab = fem::int0;
  arr<float, 2> estif(dimension(18, 18), fem::fill0);
  int kgasp = fem::int0;
  int igaus = fem::int0;
  float exisp = fem::float0;
  int jgaus = fem::int0;
  float etasp = fem::float0;
  arr_2d<4, 4, float> dmatx(fem::fill0);
  arr_2d<2, 9, float> deriv(fem::fill0);
  arr_1d<9, float> shape(fem::fill0);
  arr_2d<2, 9, float> cartd(fem::fill0);
  float djacb = fem::float0;
  arr_2d<2, 9, float> gpcod(fem::fill0);
  float dvolu = fem::float0;
  arr_2d<4, 18, float> bmatx(fem::fill0);
  int istr1 = fem::int0;
  arr_1d<4, float> stres(fem::fill0);
  arr_1d<4, float> devia(fem::fill0);
  float sint3 = fem::float0;
  float steff = fem::float0;
  float theta = fem::float0;
  float varj2 = fem::float0;
  float yield = fem::float0;
  arr_1d<4, float> avect(fem::fill0);
  float abeta = fem::float0;
  arr_1d<4, float> dvect(fem::fill0);
  int istre = fem::int0;
  int jstre = fem::int0;
  arr_2d<4, 18, float> dbmat(fem::fill0);
  //C********************************************************************
  //C
  //C        EVALUATES STIFFNESS MATRAX FOR 2-D PLANE STRESS/STRAIN
  //C        AND AXISYMMETRIC ELEMENTS
  //C
  //C********************************************************************
  twopi = 6.283185307179586f;
  cmn.io.rewind(1);
  kgaus = 0;
  //C
  //C*** LOOP OVER EACH ELMEMENT
  //C
  nstr1 = 4;
  FEM_DO_SAFE(ielem, 1, nelem) {
    lprop = matno(ielem);
    //C
    //C*** EVALUATE THE COORDINATES OF THE ELEMENT NODAL POINTS
    //C
    FEM_DO_SAFE(inode, 1, nnode) {
      lnode = fem::iabs(lnods(ielem, inode));
      iposn = (lnode - 1) * 2;
      FEM_DO_SAFE(idime, 1, 2) {
        iposn++;
        elcod(idime, inode) = coord(lnode, idime);
      }
    }
    thick = props(lprop, 3);
    //C
    //C*** INITIALIZE THE STIFFNESS MATRIX 171=NEVAB*(NEVAB+1)/2
    //C
    FEM_DO_SAFE(ievab, 1, nevab) {
      FEM_DO_SAFE(jevab, 1, nevab) {
        estif(ievab, jevab) = 0.0f;
      }
    }
    kgasp = 0;
    //C
    //C*** ENTER LOOP FOR AREA NUMERICAL INTEGRATION
    //C
    FEM_DO_SAFE(igaus, 1, ngaus) {
      exisp = posgp(igaus);
      FEM_DO_SAFE(jgaus, 1, ngaus) {
        etasp = posgp(jgaus);
        kgasp++;
        kgaus++;
        //C
        //C***  EVALUATE THE D-MATRIX
        //C
        modps(dmatx, lprop, mmats, ntype, props);
        //C
        //C***  EVALUATE THE SHAPE FUNCTIONS, ELEMENTAL VOLUME, ETC.
        //C
        sfr2(deriv, etasp, exisp, nnode, shape);
        jacob2(cmn, cartd, deriv, djacb, elcod, gpcod, ielem, kgasp,
          nnode, shape);
        dvolu = djacb * weigp(igaus) * weigp(jgaus);
        if (ntype == 3) {
          dvolu = dvolu * twopi * gpcod(1, kgasp);
        }
        if (ntype == 1) {
          dvolu = dvolu * thick;
        }
        //C
        //C*** EVALUATE THE B AND DB MATRICES
        //C
        bmatps(bmatx, cartd, nnode, shape, gpcod, ntype, kgasp);
        if (iincs == 1) {
          goto statement_80;
        }
        if (epstn(kgaus) == 0.0f) {
          goto statement_80;
        }
        FEM_DO_SAFE(istr1, 1, nstr1) {
          stres(istr1) = strsg(istr1, kgaus);
        }
        invar(devia, lprop, mmats, ncrit, props, sint3, steff, stres,
          theta, varj2, yield);
        yieldf(avect, devia, lprop, mmats, ncrit, nstr1, props,
          sint3, steff, theta, varj2);
        flowpl(avect, abeta, dvect, ntype, props, lprop, nstr1, mmats);
        FEM_DO_SAFE(istre, 1, nstre) {
          FEM_DO_SAFE(jstre, 1, nstre) {
            dmatx(istre, jstre) = dmatx(istre, jstre) - abeta * dvect(
              istre) * dvect(jstre);
          }
        }
        statement_80:
        dbe(bmatx, dbmat, dmatx, mevab, nevab, nstre, nstr1);
        //C
        //C*** CALCULATE THE ELEMENT STIFFNESSES
        //C
        FEM_DO_SAFE(ievab, 1, nevab) {
          FEM_DO_SAFE(jevab, 1, nevab) {
            FEM_DO_SAFE(istre, 1, nstre) {
              estif(ievab, jevab) += bmatx(istre, ievab) * dbmat(istre,
                jevab) * dvolu;
            }
          }
        }
      }
    }
    //C
    //C***  CONSTRUCT THE LOWER TRIANGLE OF THE STIFFNESS MATRIX
    //C
    FEM_DO_SAFE(ievab, 1, nevab) {
      FEM_DO_SAFE(jevab, 1, nevab) {
        estif(jevab, ievab) = estif(ievab, jevab);
      }
    }
    //C
    //C***  STORE THE STIFFNESS MATRIX, STRESS MATRIX AND SAMPLING POINT
    //C
    write(1, fem::unformatted), estif;
  }
}

//C$DEBUG
//C$LARGE
//C(Page 238)
void
zero(
  arr_ref<float, 2> eload,
  int const& melem,
  int const& mevab,
  int const& /* mpoin */,
  int const& mtotg,
  int const& mtotv,
  int const& ndofn,
  int const& nelem,
  int const& nevab,
  int const& /* ngaus */,
  int const& nstr1,
  int const& ntotg,
  arr_ref<float> epstn,
  arr_ref<float> effst,
  int const& ntotv,
  int const& nvfix,
  arr_ref<float, 2> strsg,
  arr_ref<float> tdisp,
  float& tfact,
  arr_ref<float, 2> tload,
  arr_ref<float, 2> treac,
  int const& mvfix)
{
  eload(dimension(melem, mevab));
  epstn(dimension(mtotg));
  effst(dimension(mtotg));
  strsg(dimension(4, mtotg));
  tdisp(dimension(mtotv));
  tload(dimension(melem, mevab));
  treac(dimension(mvfix, 2));
  //C*********************************************************************
  //C
  //C*** THIS SUBROUTINE INITIALISES VARIOUS ARRAYS TO ZERO
  //C
  //C*********************************************************************
  tfact = 0.0f;
  int ielem = fem::int0;
  int ievab = fem::int0;
  FEM_DO_SAFE(ielem, 1, nelem) {
    FEM_DO_SAFE(ievab, 1, nevab) {
      eload(ielem, ievab) = 0.0f;
      tload(ielem, ievab) = 0.0f;
    }
  }
  int itotv = fem::int0;
  FEM_DO_SAFE(itotv, 1, ntotv) {
    tdisp(itotv) = 0.0f;
  }
  int ivfix = fem::int0;
  int idofn = fem::int0;
  FEM_DO_SAFE(ivfix, 1, nvfix) {
    FEM_DO_SAFE(idofn, 1, ndofn) {
      treac(ivfix, idofn) = 0.0f;
    }
  }
  int itotg = fem::int0;
  int istr1 = fem::int0;
  FEM_DO_SAFE(itotg, 1, ntotg) {
    epstn(itotg) = 0.0f;
    effst(itotg) = 0.0f;
    FEM_DO_SAFE(istr1, 1, nstr1) {
      strsg(istr1, itotg) = 0.0f;
    }
  }
}

//C$DEBUG
//C$LARGE
//C(Page 260)
void
program_masterplast(
  int argc,
  char const* argv[])
{
  common cmn(argc, argv);
  int mbufa = fem::int0;
  int melem = fem::int0;
  int mevab = fem::int0;
  int mfron = fem::int0;
  int mmats = fem::int0;
  int mpoin = fem::int0;
  int mstif = fem::int0;
  int mtotg = fem::int0;
  int mtotv = fem::int0;
  int mvfix = fem::int0;
  int ndofn = fem::int0;
  int nprop = fem::int0;
  int nstre = fem::int0;
  arr<float, 2> coord(dimension(150, 2), fem::fill0);
  arr<int> iffix(dimension(300), fem::fill0);
  arr<int, 2> lnods(dimension(40, 9), fem::fill0);
  arr_1d<40, int> matno(fem::fill0);
  int nalgo = fem::int0;
  int ncrit = fem::int0;
  arr_1d<40, int> ndfro(fem::fill0);
  int nelem = fem::int0;
  int nevab = fem::int0;
  int ngaus = fem::int0;
  int ngau2 = fem::int0;
  int nincs = fem::int0;
  int nmats = fem::int0;
  int nnode = fem::int0;
  arr_1d<30, int> nofix(fem::fill0);
  int npoin = fem::int0;
  int nstr1 = fem::int0;
  int ntotg = fem::int0;
  int ntotv = fem::int0;
  int ntype = fem::int0;
  int nvfix = fem::int0;
  arr_1d<4, float> posgp(fem::fill0);
  arr_2d<30, 2, float> presc(fem::fill0);
  arr_2d<5, 7, float> props(fem::fill0);
  arr_1d<4, float> weigp(fem::fill0);
  arr<float, 2> rload(dimension(40, 18), fem::fill0);
  arr<float, 2> eload(dimension(40, 18), fem::fill0);
  arr<float> epstn(dimension(360), fem::fill0);
  arr<float> effst(dimension(360), fem::fill0);
  arr<float, 2> strsg(dimension(4, 360), fem::fill0);
  arr<float> tdisp(dimension(300), fem::fill0);
  float tfact = fem::float0;
  arr<float, 2> tload(dimension(40, 18), fem::fill0);
  arr_2d<30, 2, float> treac(fem::fill0);
  int iincs = fem::int0;
  arr<float> fixed(dimension(300), fem::fill0);
  int miter = fem::int0;
  arr_1d<2, int> noutp(fem::fill0);
  float toler = fem::float0;
  int iiter = fem::int0;
  int kresl = fem::int0;
  arr<float> asdis(dimension(300), fem::fill0);
  arr_1d<10, float> eqrhs(fem::fill0);
  arr<float, 2> equat(dimension(80, 10), fem::fill0);
  arr<float, 2> estif(dimension(18, 18), fem::fill0);
  arr_1d<80, float> gload(fem::fill0);
  arr<float> gstif(dimension(3240), fem::fill0);
  arr_1d<18, int> locel(fem::fill0);
  arr_1d<80, int> nacva(fem::fill0);
  arr_1d<10, int> namev(fem::fill0);
  arr_1d<18, int> ndest(fem::fill0);
  arr_1d<10, int> npivo(fem::fill0);
  arr_1d<80, float> vecrv(fem::fill0);
  float facto = fem::float0;
  int lprop = fem::int0;
  int nchek = fem::int0;
  float pvalu = fem::float0;
  arr<float> stfor(dimension(300), fem::fill0);
  arr<float> tofor(dimension(300), fem::fill0);
  //C****************************************************************
  //C       PROGRAM FOR THE ELASTO-PLASTIC ANALYSIS OF PLANE STRESS,
  //C       PLANE STRAIN END AXISYMMETRIC SOLIDS
  //C****************************************************************
  //C***
  cmn.io.open(5, "C:\\Users\\zhangzhg\\Desktop\\FEM\\Debug\\data\\input.dat")
    .status("unknown");
  cmn.io.open(6, "C:\\Users\\zhangzhg\\Desktop\\FEM\\Debug\\data\\output.dat")
    .status("unknown");
  //C***
  //C
  //C*** PRESET VARIABLES ASSOCIATED WITH DYNAMIC DIMENSIONING
  //C
  dimen(mbufa, melem, mevab, mfron, mmats, mpoin, mstif, mtotg,
    mtotv, mvfix, ndofn, nprop, nstre);
  //C
  //C*** CALL THE SUBRIUTINE WHICH READS MOST OF PROBLEM DATA
  //C
  input(cmn, coord, iffix, lnods, matno, melem, mevab, mfron, mmats,
    mpoin, mtotv, mvfix, nalgo, ncrit, ndfro, ndofn, nelem, nevab,
    ngaus, ngau2, nincs, nmats, nnode, nofix, npoin, nprop, nstre,
    nstr1, ntotg, ntotv, ntype, nvfix, posgp, presc, props, weigp);
  //C
  //C*** CALL THE SUBROUTINE WHICH COMPUTES THE CONSISTENT LOAD VECTORS
  //C    FOR EACH ELEMENT AFTER READING THE RELEVENT INPUT DATA
  //C
  loadps(cmn, coord, lnods, matno, melem, mmats, mpoin, nelem, nevab,
    ngaus, nnode, npoin, nstre, ntype, posgp, props, rload, weigp,
    ndofn);
  //C
  //C*** INITIALISE CERTAIN ARRAYS
  //C
  zero(eload, melem, mevab, mpoin, mtotg, mtotv, ndofn, nelem, nevab,
    ngaus, nstr1, ntotg, epstn, effst, ntotv, nvfix, strsg, tdisp,
    tfact, tload, treac, mvfix);
  //C
  //C*** LOOP OVER EACH INCREMENT
  //C
  FEM_DO_SAFE(iincs, 1, nincs) {
    //C
    //C*** READ DATA FOR CURRENT INCREMENT
    //C
    increm(cmn, eload, fixed, iincs, melem, mevab, miter, mtotv,
      mvfix, ndofn, nelem, nevab, noutp, nofix, ntotv, nvfix, presc,
      rload, tfact, tload, toler);
    //C
    //C*** LOOP OVER EACH ITERATION
    //C
    FEM_DO_SAFE(iiter, 1, miter) {
      //C
      //C*** CALL ROUTINE WHICH SELECTS SOLUTION ALORITHM VARIABLE KRESL
      //C
      algor(fixed, iincs, iiter, kresl, mtotv, nalgo, ntotv);
      //C
      //C*** CHECK WHETHER A NEW EVALUATION OF THE STIFFNESS MATRIX IS REQUIRED
      //C
      if (kresl == 1) {
        stiffp(cmn, coord, epstn, iincs, lnods, matno, mevab, mmats,
          mpoin, mtotv, nelem, nevab, ngaus, nnode, nstre, nstr1,
          posgp, props, weigp, melem, mtotg, strsg, ntype, ncrit);
      }
      //C
      //C*** SOLVE EQUATIONS
      //C
      front(cmn, asdis, eload, eqrhs, equat, estif, fixed, iffix,
        iincs, iiter, gload, gstif, locel, lnods, kresl, mbufa,
        melem, mevab, mfron, mstif, mtotv, mvfix, nacva, namev,
        ndest, ndofn, nelem, nevab, nnode, nofix, npivo, npoin,
        ntotv, tdisp, tload, treac, vecrv);
      //C
      //C*** CALCULATE RESIDUAL FORCES
      //C
      residu(cmn, asdis, coord, effst, eload, facto, iiter, lnods,
        lprop, matno, melem, mmats, mpoin, mtotg, mtotv, ndofn, nelem,
        nevab, ngaus, nnode, nstr1, ntype, posgp, props, nstre, ncrit,
        strsg, weigp, tdisp, epstn);
      //C
      //C*** CHECK FOR CONVERGENCE
      //C
      conver(cmn, eload, iiter, lnods, melem, mevab, mtotv, nchek,
        ndofn, nelem, nevab, nnode, ntotv, pvalu, stfor, tload,
        tofor, toler);
      //C
      //C*** OUTPUT RESULTS IF REQUIRED
      //C
      if (iiter == 1 && noutp(1) > 0) {
        output(cmn, iiter, mtotg, mtotv, mvfix, nelem, ngaus, nofix,
          noutp, npoin, nvfix, strsg, tdisp, treac, epstn, ntype,
          nchek);
      }
      //C
      //C*** IF SOLUTION HAS CONVERGED STOP ITERATING AND OUTPUT RESULTS
      //C
      if (nchek == 0) {
        goto statement_75;
      }
    }
    //C
    //C***
    //C
    if (nalgo == 2) {
      goto statement_75;
    }
    FEM_STOP(0);
    statement_75:
    output(cmn, iiter, mtotg, mtotv, mvfix, nelem, ngaus, nofix,
      noutp, npoin, nvfix, strsg, tdisp, treac, epstn, ntype, nchek);
  }
  //C***
  cmn.io.close(5);
  cmn.io.close(6);
  //C***
  FEM_STOP(0);
}

} // namespace epca

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    epca::program_masterplast);
}
