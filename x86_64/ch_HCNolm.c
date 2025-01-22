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
 
#define nrn_init _nrn_init__ch_HCNolm
#define _nrn_initial _nrn_initial__ch_HCNolm
#define nrn_cur _nrn_cur__ch_HCNolm
#define _nrn_current _nrn_current__ch_HCNolm
#define nrn_jacob _nrn_jacob__ch_HCNolm
#define nrn_state _nrn_state__ch_HCNolm
#define _net_receive _net_receive__ch_HCNolm 
#define _f_rates _f_rates__ch_HCNolm 
#define deriv deriv__ch_HCNolm 
#define rates rates__ch_HCNolm 
 
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
#define g _p[1]
#define g_columnindex 1
#define ih _p[2]
#define ih_columnindex 2
#define myi _p[3]
#define myi_columnindex 3
#define r _p[4]
#define r_columnindex 4
#define eh _p[5]
#define eh_columnindex 5
#define Dr _p[6]
#define Dr_columnindex 6
#define _g _p[7]
#define _g_columnindex 7
#define _ion_eh	*_ppvar[0]._pval
#define _ion_ih	*_ppvar[1]._pval
#define _ion_dihdv	*_ppvar[2]._pval
 
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
 static void _hoc_efun(void);
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
 "setdata_ch_HCNolm", _hoc_setdata,
 "efun_ch_HCNolm", _hoc_efun,
 "rates_ch_HCNolm", _hoc_rates,
 0, 0
};
#define efun efun_ch_HCNolm
 extern double efun( double );
 /* declare global and static user variables */
#define rexp rexp_ch_HCNolm
 double rexp = 0;
#define rinf rinf_ch_HCNolm
 double rinf = 0;
#define tau_r tau_r_ch_HCNolm
 double tau_r = 0;
#define usetable usetable_ch_HCNolm
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_ch_HCNolm", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_ch_HCNolm", "mho/cm2",
 "g_ch_HCNolm", "mho/cm2",
 "ih_ch_HCNolm", "mA/cm2",
 "myi_ch_HCNolm", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double r0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "rinf_ch_HCNolm", &rinf_ch_HCNolm,
 "rexp_ch_HCNolm", &rexp_ch_HCNolm,
 "tau_r_ch_HCNolm", &tau_r_ch_HCNolm,
 "usetable_ch_HCNolm", &usetable_ch_HCNolm,
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
"ch_HCNolm",
 "gmax_ch_HCNolm",
 "g_ch_HCNolm",
 0,
 "ih_ch_HCNolm",
 "myi_ch_HCNolm",
 0,
 "r_ch_HCNolm",
 0,
 0};
 static Symbol* _h_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	gmax = 0.001385;
 	g = 0;
 	_prop->param = _p;
 	_prop->param_size = 8;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_h_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* eh */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ih */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dihdv */
 
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

 void _ch_HCNolm_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("h", 1.0);
 	_h_sym = hoc_lookup("h_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 8, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "h_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ch_HCNolm /home/dimitrios/TFM/CA1model/mechanisms/ch_HCNolm.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_rinf;
 static double *_t_tau_r;
 static double *_t_rexp;
static int _reset;
static char *modelname = "Hyperpolarization-activated, CN-gated channel (voltage dependent, for O-LM cells)";

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
 static int _slist2[1]; static double _dlist2[1];
 static double _savstate1[1], *_temp1 = _savstate1;
 static int _slist1[1], _dlist1[1];
 static int deriv(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dr = ( rinf - r ) / tau_r ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dr = Dr  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_r )) ;
  return 0;
}
 /*END CVODE*/
 
static int deriv () {_reset=0;
 { static int _recurse = 0;
 int _counte = -1;
 if (!_recurse) {
 _recurse = 1;
 {int _id; for(_id=0; _id < 1; _id++) { _savstate1[_id] = _p[_slist1[_id]];}}
 error = newton(1,_slist2, _p, deriv, _dlist2);
 _recurse = 0; if(error) {abort_run(error);}}
 {
   rates ( _threadargscomma_ v ) ;
   Dr = ( rinf - r ) / tau_r ;
   {int _id; for(_id=0; _id < 1; _id++) {
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
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 200.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/300.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 301; _x += _dx, _i++) {
    _f_rates(_x);
    _t_rinf[_i] = rinf;
    _t_tau_r[_i] = tau_r;
    _t_rexp[_i] = rexp;
   }
   _sav_dt = dt;
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
  rinf = _xi;
  tau_r = _xi;
  rexp = _xi;
  return;
 }
 if (_xi <= 0.) {
 rinf = _t_rinf[0];
 tau_r = _t_tau_r[0];
 rexp = _t_rexp[0];
 return; }
 if (_xi >= 300.) {
 rinf = _t_rinf[300];
 tau_r = _t_tau_r[300];
 rexp = _t_rexp[300];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 rinf = _t_rinf[_i] + _theta*(_t_rinf[_i+1] - _t_rinf[_i]);
 tau_r = _t_tau_r[_i] + _theta*(_t_tau_r[_i+1] - _t_tau_r[_i]);
 rexp = _t_rexp[_i] + _theta*(_t_rexp[_i+1] - _t_rexp[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   rinf = 1.0 / ( 1.0 + exp ( ( _lv + 84.1 ) / 10.2 ) ) ;
   tau_r = 100.0 + 1.0 / ( exp ( - 17.9 - 0.116 * _lv ) + exp ( - 1.84 + 0.09 * _lv ) ) ;
   rexp = 1.0 - exp ( - dt / ( tau_r ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
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
  eh = _ion_eh;
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
  eh = _ion_eh;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_h_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_h_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_h_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  r = r0;
 {
   rates ( _threadargscomma_ v ) ;
   r = rinf ;
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
  eh = _ion_eh;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * r ;
   ih = g * ( v - eh ) ;
   myi = ih ;
   }
 _current += ih;

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
  eh = _ion_eh;
 _g = _nrn_current(_v + .001);
 	{ double _dih;
  _dih = ih;
 _rhs = _nrn_current(_v);
  _ion_dihdv += (_dih - ih)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ih += ih ;
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
  eh = _ion_eh;
 { error = _deriv1_advance = 1;
 derivimplicit(_ninits, 1, _slist1, _dlist1, _p, &t, dt, deriv, &_temp1);
_deriv1_advance = 0;
 if(error){fprintf(stderr,"at line 91 in file ch_HCNolm.mod:\n	SOLVE deriv METHOD derivimplicit\n"); nrn_complain(_p); abort_run(error);}
    if (secondorder) {
    int _i;
    for (_i = 0; _i < 1; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = r_columnindex;  _dlist1[0] = Dr_columnindex;
 _slist2[0] = r_columnindex;
   _t_rinf = makevector(301*sizeof(double));
   _t_tau_r = makevector(301*sizeof(double));
   _t_rexp = makevector(301*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/dimitrios/TFM/CA1model/mechanisms/ch_HCNolm.mod";
static const char* nmodl_file_text = 
  "TITLE Hyperpolarization-activated, CN-gated channel (voltage dependent, for O-LM cells)\n"
  "\n"
  "COMMENT\n"
  "Hyperpolarization-activated, CN-gated channel (voltage dependent, for O-LM cells)\n"
  "\n"
  "Ions: non-specific\n"
  "\n"
  "Style: quasi-ohmic\n"
  "\n"
  "From: 1.	Maccaferri, G. and McBain, C.J. The hyperpolarization-activated current\n"
  "	(Ih) and its contribution to pacemaker activity in rat CA1 hippocampal\n"
  "	stratum oriens-alveus interneurons, J. Physiol. 497.1:119-130,\n"
  "	1996.\n"
  "\n"
  "		V1/2 = -84.1 mV\n"
  "		   k = 10.2\n"
  "		reversal potential = -32.9 +/- 1.1 mV\n"
  "\n"
  "at -70 mV, currents were fitted by a single exponetial of t = 2.8+/- 0.76 s\n"
  "at -120 mV, two exponentials were required, t1 = 186.3+/-33.6 ms \n"
  "t2 = 1.04+/-0.16 s\n"
  "\n"
  "\n"
  "2.	Maccaferri, G. et al. Properties of the\n"
  "	Hyperpoarization-activated current in rat hippocampal CA1 Pyramidal\n"
  "	cells. J. Neurophysiol. Vol. 69 No. 6:2129-2136, 1993.\n"
  "\n"
  "		V1/2 = -97.9 mV\n"
  "		   k = 13.4\n"
  "		reversal potential = -18.3 mV\n"
  "\n"
  "3.	Pape, H.C.  Queer current and pacemaker: The\n"
  "	hyperpolarization-activated cation current in neurons, Annu. Rev. \n"
  "	Physiol. 58:299-327, 1996.\n"
  "\n"
  "		single channel conductance is around 1 pS\n"
  "		average channel density is below 0.5 um-2\n"
  "		0.5 pS/um2 = 0.00005 mho/cm2 = 0.05 umho/cm2		\n"
  "4.	Magee, J.C. Dendritic Hyperpolarization-Activated Currents Modify\n"
  "	the Integrative Properties of Hippocampal CA1 Pyramidal Neurons, J.\n"
  "	Neurosci., 18(19):7613-7624, 1998\n"
  "\n"
  "Deals with Ih in CA1 pyramidal cells.  Finds that conductance density\n"
  "increases with distance from the soma.\n"
  "\n"
  "soma g = 0.0013846 mho/cm2\n"
  "dendrite g (300-350 um away) = 0.0125 mho/cm2\n"
  "see Table 1 in th paper\n"
  "\n"
  "Updates:\n"
  "2014 December (Marianne Bezaire): documented\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  " UNITS {\n"
  "        (mA) = (milliamp)\n"
  "        (mV) = (millivolt)\n"
  "}\n"
  " \n"
  "NEURON {\n"
  "        SUFFIX ch_HCNolm\n"
  "        USEION h READ eh WRITE ih VALENCE 1\n"
  "        RANGE gmax,ih, g\n"
  "        GLOBAL rinf, rexp, tau_r\n"
  "        RANGE myi\n"
  "}\n"
  " \n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  " \n"
  "PARAMETER {\n"
  "        v (mV)\n"
  "        dt (ms)\n"
  "        gmax = 0.001385 (mho/cm2)\n"
  "		g (mho/cm2)\n"
  "        eh = -32.9 (mV)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "        r\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "        ih (mA/cm2)\n"
  "	rinf rexp\n"
  "	tau_r\n"
  "	myi (mA/cm2)\n"
  "}\n"
  " \n"
  "BREAKPOINT {\n"
  "	SOLVE deriv METHOD derivimplicit\n"
  "	g = gmax*r\n"
  "	ih = g*(v - eh)\n"
  "	myi = ih\n"
  "}\n"
  " \n"
  "INITIAL {\n"
  "	rates(v)\n"
  "	r = rinf\n"
  "}\n"
  "\n"
  "DERIVATIVE deriv { :Computes state variable h at current v and dt.\n"
  "	rates(v)\n"
  "	r' = (rinf - r)/tau_r\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "        TABLE rinf, tau_r, rexp DEPEND dt FROM -200\n"
  "TO 100 WITH 300\n"
  "	rinf = 1/(1 + exp((v+84.1)/10.2))\n"
  "	tau_r = 100 + 1/(exp(-17.9-0.116*v)+exp(-1.84+0.09*v))\n"
  "	rexp = 1 - exp(-dt/(tau_r))\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "  if (fabs(z) < 1e-4) {\n"
  "    efun = 1 - z/2\n"
  "  } else {\n"
  "    efun = z/(exp(z) - 1)\n"
  "  }\n"
  "}\n"
  "\n"
  " \n"
  "UNITSON\n"
  "\n"
  ;
#endif
