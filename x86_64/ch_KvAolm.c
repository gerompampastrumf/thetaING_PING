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
 
#define nrn_init _nrn_init__ch_KvAolm
#define _nrn_initial _nrn_initial__ch_KvAolm
#define nrn_cur _nrn_cur__ch_KvAolm
#define _nrn_current _nrn_current__ch_KvAolm
#define nrn_jacob _nrn_jacob__ch_KvAolm
#define nrn_state _nrn_state__ch_KvAolm
#define _net_receive _net_receive__ch_KvAolm 
#define _f_rates _f_rates__ch_KvAolm 
#define deriv deriv__ch_KvAolm 
#define rates rates__ch_KvAolm 
 
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
#define ik _p[1]
#define ik_columnindex 1
#define myi _p[2]
#define myi_columnindex 2
#define g _p[3]
#define g_columnindex 3
#define a _p[4]
#define a_columnindex 4
#define b _p[5]
#define b_columnindex 5
#define ek _p[6]
#define ek_columnindex 6
#define Da _p[7]
#define Da_columnindex 7
#define Db _p[8]
#define Db_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 /* declaration of user functions */
 static void _hoc_rates(void);
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
 "setdata_ch_KvAolm", _hoc_setdata,
 "rates_ch_KvAolm", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define aexp aexp_ch_KvAolm
 double aexp = 0;
#define ainf ainf_ch_KvAolm
 double ainf = 0;
#define bexp bexp_ch_KvAolm
 double bexp = 0;
#define binf binf_ch_KvAolm
 double binf = 0;
#define p p_ch_KvAolm
 double p = 5;
#define tau_a tau_a_ch_KvAolm
 double tau_a = 5;
#define tau_b tau_b_ch_KvAolm
 double tau_b = 0;
#define usetable usetable_ch_KvAolm
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_ch_KvAolm", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "p_ch_KvAolm", "degC",
 "tau_a_ch_KvAolm", "ms",
 "gmax_ch_KvAolm", "mho/cm2",
 "ik_ch_KvAolm", "mA/cm2",
 "myi_ch_KvAolm", "mA/cm2",
 "g_ch_KvAolm", "mho/cm2",
 0,0
};
 static double a0 = 0;
 static double b0 = 0;
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "p_ch_KvAolm", &p_ch_KvAolm,
 "tau_a_ch_KvAolm", &tau_a_ch_KvAolm,
 "ainf_ch_KvAolm", &ainf_ch_KvAolm,
 "binf_ch_KvAolm", &binf_ch_KvAolm,
 "aexp_ch_KvAolm", &aexp_ch_KvAolm,
 "bexp_ch_KvAolm", &bexp_ch_KvAolm,
 "tau_b_ch_KvAolm", &tau_b_ch_KvAolm,
 "usetable_ch_KvAolm", &usetable_ch_KvAolm,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ch_KvAolm",
 "gmax_ch_KvAolm",
 0,
 "ik_ch_KvAolm",
 "myi_ch_KvAolm",
 "g_ch_KvAolm",
 0,
 "a_ch_KvAolm",
 "b_ch_KvAolm",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gmax = 0.0165;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ch_KvAolm_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ch_KvAolm /home/dimitrios/TFM/CA1model/mechanisms/ch_KvAolm.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_ainf;
 static double *_t_aexp;
 static double *_t_binf;
 static double *_t_bexp;
 static double *_t_tau_a;
 static double *_t_tau_b;
static int _reset;
static char *modelname = "A-type potassium channel (voltage dependent, for O-LM cell)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 static int _deriv1_advance = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist2[2]; static double _dlist2[2];
 static double _savstate1[2], *_temp1 = _savstate1;
 static int _slist1[2], _dlist1[2];
 static int deriv(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   Da = ( ainf - a ) / ( tau_a ) ;
   Db = ( binf - b ) / ( tau_b ) ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 Da = Da  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( tau_a ) )) ;
 Db = Db  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ( tau_b ) )) ;
  return 0;
}
 /*END CVODE*/
 
static int deriv () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 2; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(2,_slist2, _p, deriv, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   Da = ( ainf - a ) / ( tau_a ) ;
   Db = ( binf - b ) / ( tau_b ) ;
   {int _id; for(_id=0; _id < 2; _id++) {
if (_deriv1_advance) {
 _dlist2[++_counte] = _p[_dlist1[_id]] - (_p[_slist1[_id]] - _savstate1[_id])/dt;
 }else{
_dlist2[++_counte] = _p[_slist1[_id]] - _savstate1[_id];}}}
 } }
 return _reset;}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_p;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_p != p) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 200.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/300.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 301; _x += _dx, _i++) {
    _f_rates(_x);
    _t_ainf[_i] = ainf;
    _t_aexp[_i] = aexp;
    _t_binf[_i] = binf;
    _t_bexp[_i] = bexp;
    _t_tau_a[_i] = tau_a;
    _t_tau_b[_i] = tau_b;
   }
   _sav_dt = dt;
   _sav_p = p;
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  ainf = _xi;
  aexp = _xi;
  binf = _xi;
  bexp = _xi;
  tau_a = _xi;
  tau_b = _xi;
  return;
 }
 if (_xi <= 0.) {
 ainf = _t_ainf[0];
 aexp = _t_aexp[0];
 binf = _t_binf[0];
 bexp = _t_bexp[0];
 tau_a = _t_tau_a[0];
 tau_b = _t_tau_b[0];
 return; }
 if (_xi >= 300.) {
 ainf = _t_ainf[300];
 aexp = _t_aexp[300];
 binf = _t_binf[300];
 bexp = _t_bexp[300];
 tau_a = _t_tau_a[300];
 tau_b = _t_tau_b[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ainf = _t_ainf[_i] + _theta*(_t_ainf[_i+1] - _t_ainf[_i]);
 aexp = _t_aexp[_i] + _theta*(_t_aexp[_i+1] - _t_aexp[_i]);
 binf = _t_binf[_i] + _theta*(_t_binf[_i+1] - _t_binf[_i]);
 bexp = _t_bexp[_i] + _theta*(_t_bexp[_i+1] - _t_bexp[_i]);
 tau_a = _t_tau_a[_i] + _theta*(_t_tau_a[_i+1] - _t_tau_a[_i]);
 tau_b = _t_tau_b[_i] + _theta*(_t_tau_b[_i+1] - _t_tau_b[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   double _lalpha_b , _lbeta_b ;
 _lalpha_b = 0.000009 / exp ( ( _lv - 26.0 ) / 18.5 ) ;
   _lbeta_b = 0.014 / ( exp ( ( _lv + 70.0 ) / ( - 11.0 ) ) + 0.2 ) ;
   ainf = 1.0 / ( 1.0 + exp ( - ( _lv + 14.0 ) / 16.6 ) ) ;
   aexp = 1.0 - exp ( - dt / ( tau_a ) ) ;
   tau_b = 1.0 / ( _lalpha_b + _lbeta_b ) ;
   binf = 1.0 / ( 1.0 + exp ( ( _lv + 71.0 ) / 7.3 ) ) ;
   bexp = 1.0 - exp ( - dt / ( tau_b ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  a = a0;
  b = b0;
 {
   rates ( _threadargscomma_ v ) ;
   a = ainf ;
   b = binf ;
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
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * a * b ;
   ik = g * ( v - ek ) ;
   myi = ik ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
  ek = _ion_ek;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 2, _slist1, _dlist1, _p, &t, dt, deriv, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 84 in file ch_KvAolm.mod:\n        SOLVE deriv METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 2; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = a_columnindex;  _dlist1[0] = Da_columnindex;
 _slist1[1] = b_columnindex;  _dlist1[1] = Db_columnindex;
 _slist2[0] = a_columnindex;
 _slist2[1] = b_columnindex;
   _t_ainf = makevector(301*sizeof(double));
   _t_aexp = makevector(301*sizeof(double));
   _t_binf = makevector(301*sizeof(double));
   _t_bexp = makevector(301*sizeof(double));
   _t_tau_a = makevector(301*sizeof(double));
   _t_tau_b = makevector(301*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/dimitrios/TFM/CA1model/mechanisms/ch_KvAolm.mod";
static const char* nmodl_file_text = 
  "TITLE A-type potassium channel (voltage dependent, for O-LM cell)\n"
  "\n"
  "COMMENT\n"
  "A-type potassium channel (voltage dependent, for O-LM cell)\n"
  "\n"
  "Ions: k\n"
  "\n"
  "Style: quasi-ohmic\n"
  "\n"
  "From: 1.	Zhang, L. and McBain, J. Voltage-gated potassium currents in\n"
  "	stratum oriens-alveus inhibitory neurons of the rat CA1\n"
  "	hippocampus, J. Physiol. 488.3:647-660, 1995.\n"
  "\n"
  "		Activation V1/2 = -14 mV\n"
  "		slope = 16.6\n"
  "		activation t = 5 ms\n"
  "		Inactivation V1/2 = -71 mV\n"
  "		slope = 7.3\n"
  "		inactivation t = 15 ms\n"
  "		recovery from inactivation = 142 ms\n"
  "\n"
  "2.	Martina, M. et al. Functional and Molecular Differences between\n"
  "	Voltage-gated K+ channels of fast-spiking interneurons and pyramidal\n"
  "	neurons of rat hippocampus, J. Neurosci. 18(20):8111-8125, 1998.	\n"
  "	(only the gmax is from this paper)\n"
  "\n"
  "		gmax = 0.0175 mho/cm2\n"
  "		Activation V1/2 = -6.2 +/- 3.3 mV\n"
  "		slope = 23.0 +/- 0.7 mV\n"
  "		Inactivation V1/2 = -75.5 +/- 2.5 mV\n"
  "		slope = 8.5 +/- 0.8 mV\n"
  "		recovery from inactivation t = 165 +/- 49 ms  \n"
  "\n"
  "3.	Warman, E.N. et al.  Reconstruction of Hippocampal CA1 pyramidal\n"
  "	cell electrophysiology by computer simulation, J. Neurophysiol.\n"
  "	71(6):2033-2045, 1994.\n"
  "\n"
  "		gmax = 0.01 mho/cm2\n"
  "		(number taken from the work by Numann et al. in guinea pig\n"
  "		CA1 neurons)\n"
  "\n"
  "Updates:\n"
  "2014 December (Marianne Bezaire): documented\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "UNITS {\n"
  "        (mA) = (milliamp)\n"
  "        (mV) = (millivolt)\n"
  "}\n"
  " \n"
  "NEURON {\n"
  "        SUFFIX ch_KvAolm\n"
  "        USEION k READ ek WRITE ik\n"
  "        RANGE gmax,ik\n"
  "        GLOBAL ainf, binf, aexp, bexp, tau_b\n"
  "        RANGE myi, g\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV)\n"
  "        p = 5 (degC)\n"
  "        dt (ms)\n"
  "        gmax = 0.0165 (mho/cm2)	:from Martina et al.\n"
  "        ek = -90 (mV)\n"
  "	tau_a = 5 (ms)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "        a b\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "        ik (mA/cm2)\n"
  "	ainf binf aexp bexp\n"
  "	tau_b\n"
  "	myi (mA/cm2)\n"
  "	g (mho/cm2)\n"
  "}\n"
  " \n"
  "BREAKPOINT {\n"
  "        SOLVE deriv METHOD derivimplicit\n"
  "		g = gmax*a*b\n"
  "        ik = g*(v - ek)\n"
  "	myi = ik\n"
  "}\n"
  " \n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	a = ainf\n"
  "	b = binf\n"
  "}\n"
  "\n"
  "DERIVATIVE deriv {  :Computes state variables m, h, and n rates(v)      \n"
  "		: at the current v and dt.\n"
  "        a' = (ainf - a)/(tau_a)\n"
  "        b' = (binf - b)/(tau_b)\n"
  "}\n"
  " \n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "        LOCAL alpha_b, beta_b\n"
  "	TABLE ainf, aexp, binf, bexp, tau_a, tau_b  DEPEND dt, p FROM -200\n"
  "TO 100 WITH 300\n"
  "	alpha_b = 0.000009/exp((v-26)/18.5)\n"
  "	beta_b = 0.014/(exp((v+70)/(-11))+0.2)\n"
  "        ainf = 1/(1 + exp(-(v + 14)/16.6))\n"
  "        aexp = 1 - exp(-dt/(tau_a))\n"
  "	tau_b = 1/(alpha_b + beta_b)\n"
  "        binf = 1/(1 + exp((v + 71)/7.3))\n"
  "        bexp = 1 - exp(-dt/(tau_b))\n"
  "}\n"
  " \n"
  "UNITSON\n"
  "\n"
  ;
#endif
