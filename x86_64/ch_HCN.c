/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__ch_HCN
#define _nrn_initial _nrn_initial__ch_HCN
#define nrn_cur _nrn_cur__ch_HCN
#define _nrn_current _nrn_current__ch_HCN
#define nrn_jacob _nrn_jacob__ch_HCN
#define nrn_state _nrn_state__ch_HCN
#define _net_receive _net_receive__ch_HCN 
#define _f_trates _f_trates__ch_HCN 
#define states states__ch_HCN 
#define trates trates__ch_HCN 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _p[0]
#define gmax_columnindex 0
#define e _p[1]
#define e_columnindex 1
#define g _p[2]
#define g_columnindex 2
#define i _p[3]
#define i_columnindex 3
#define hinf _p[4]
#define hinf_columnindex 4
#define slow_tau _p[5]
#define slow_tau_columnindex 5
#define myi _p[6]
#define myi_columnindex 6
#define h _p[7]
#define h_columnindex 7
#define Dh _p[8]
#define Dh_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_trates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ch_HCN", _hoc_setdata,
 "trates_ch_HCN", _hoc_trates,
 0, 0
};
 /* declare global and static user variables */
#define usetable usetable_ch_HCN
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_ch_HCN", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_ch_HCN", "mho/cm2",
 "e_ch_HCN", "mV",
 "g_ch_HCN", "mho/cm2",
 "i_ch_HCN", "mA/cm2",
 "slow_tau_ch_HCN", "ms",
 "myi_ch_HCN", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_ch_HCN", &usetable_ch_HCN,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ch_HCN",
 "gmax_ch_HCN",
 "e_ch_HCN",
 0,
 "g_ch_HCN",
 "i_ch_HCN",
 "hinf_ch_HCN",
 "slow_tau_ch_HCN",
 "myi_ch_HCN",
 0,
 "h_ch_HCN",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	e = 0;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ch_HCN_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ch_HCN /home/dimitrios/TFM/CA1model/mechanisms/ch_HCN.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double _zq10 ;
 static double *_t_hinf;
 static double *_t_slow_tau;
static int _reset;
static char *modelname = "Hyperpolarization-activated, CN-gated channel (voltage dependent)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int trates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_trates(double);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*VERBATIM*/
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
   Dh = ( hinf - h ) / slow_tau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v ) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / slow_tau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / slow_tau)))*(- ( ( ( hinf ) ) / slow_tau ) / ( ( ( ( - 1.0 ) ) ) / slow_tau ) - h) ;
   }
  return 0;
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 120.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_trates)/220.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 221; _x += _dx, _i++) {
    _f_trates(_x);
    _t_hinf[_i] = hinf;
    _t_slow_tau[_i] = slow_tau;
   }
   _sav_celsius = celsius;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  hinf = _xi;
  slow_tau = _xi;
  return;
 }
 if (_xi <= 0.) {
 hinf = _t_hinf[0];
 slow_tau = _t_slow_tau[0];
 return; }
 if (_xi >= 220.) {
 hinf = _t_hinf[220];
 slow_tau = _t_slow_tau[220];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 slow_tau = _t_slow_tau[_i] + _theta*(_t_slow_tau[_i+1] - _t_slow_tau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   _zq10 = pow( 3.0 , ( ( celsius - 34.0 ) / 10.0 ) ) ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv + 91.0 ) / 10.0 ) ) ;
   slow_tau = ( 80.0 * 1.5 + .75 * 172.7 / ( 1.0 + exp ( - ( _lv + 59.3 ) / - 0.83 ) ) ) / _zq10 ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
 {
   trates ( _threadargscomma_ v ) ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * h * h ;
   i = g * ( v - e ) ;
   myi = i ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  states();
 if(error){fprintf(stderr,"at line 72 in file ch_HCN.mod:\n		\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = h_columnindex;  _dlist1[0] = Dh_columnindex;
   _t_hinf = makevector(221*sizeof(double));
   _t_slow_tau = makevector(221*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/dimitrios/TFM/CA1model/mechanisms/ch_HCN.mod";
static const char* nmodl_file_text = 
  "TITLE Hyperpolarization-activated, CN-gated channel (voltage dependent)\n"
  "\n"
  "COMMENT\n"
  "Hyperpolarization-activated, CN-gated channel (voltage dependent)\n"
  "\n"
  "Ions: non-specific\n"
  "\n"
  "Style: quasi-ohmic\n"
  "\n"
  "From: Chen et al (2001), distal dendrite, control h channel\n"
  "(modeling by Ildiko Aradi, iaradi@uci.edu)\n"
  "\n"
  "Updates:\n"
  "2014 December (Marianne Bezaire): documented, only slow tau used, fake\n"
  "ion switched to a non-specific current\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "VERBATIM\n"
  "#include <stdlib.h> /* 	Include this library so that the following\n"
  "						(innocuous) warning does not appear:\n"
  "						 In function '_thread_cleanup':\n"
  "						 warning: incompatible implicit declaration of \n"
  "						          built-in function 'free'  */\n"
  "ENDVERBATIM\n"
  " \n"
  "UNITS {\n"
  "	(mA) =(milliamp)\n"
  "	(mV) =(millivolt)\n"
  "	(uF) = (microfarad)\n"
  "	(molar) = (1/liter)\n"
  "	(nA) = (nanoamp)\n"
  "	(mM) = (millimolar)\n"
  "	(um) = (micron)\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134	(joule/degC)\n"
  "}\n"
  " \n"
  "NEURON { \n"
  "	SUFFIX ch_HCN \n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE gmax, g, i, e\n"
  "	RANGE hinf\n"
  "	RANGE slow_tau :, fast_tau\n"
  "	RANGE myi\n"
  "}\n"
  " \n"
  " \n"
  "PARAMETER {\n"
  "	gmax  (mho/cm2)\n"
  "	e     (mV)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	h\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "	v         (mV) \n"
  "	celsius   (degC)\n"
  "	dt        (ms)    \n"
  "	g         (mho/cm2)\n"
  " 	i         (mA/cm2)\n"
  "	hinf \n"
  " 	:fast_tau (ms)\n"
  " 	slow_tau  (ms) \n"
  " 	myi       (mA/cm2)\n"
  "} \n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "		\n"
  "	g = gmax*h*h\n"
  "	i = g*(v - e)\n"
  "	myi = i\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL { : called from hoc to calculate hinf at resting potential\n"
  "	trates(v)\n"
  "	h = hinf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {	:computes h at current v and dt \n"
  "	trates(v)\n"
  "	h' = (hinf-h)/slow_tau :  + (hinf-h)/fast_tau\n"
  "}\n"
  " \n"
  "LOCAL q10\n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "	TABLE hinf, slow_tau : , fast_tau\n"
  "	DEPEND celsius \n"
  "	FROM -120 TO 100 WITH 220\n"
  "                           \n"
  "    :q10 = 3^((celsius - 6.3)/10)\n"
  "    q10 = 3^((celsius - 34)/10)\n"
  "       \n"
  "	hinf =  1 / (1 + exp( (v+91)/10 ))\n"
  "\n"
  "	:\"hyf\" FAST CONTROL Hype activation system\n"
  "	:fast_tau = (14.9 + 14.1 / (1+exp(-(v+95.2)/0.5)))/q10\n"
  "\n"
  "	:\"hys\" SLOW CONTROL Hype activation system\n"
  "	slow_tau = (80*1.5 + .75*172.7 / (1+exp(-(v+59.3)/-0.83)))/q10\n"
  "\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
