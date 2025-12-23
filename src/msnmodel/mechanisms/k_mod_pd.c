/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
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
 
#define nrn_init _nrn_init__k_mod_pd
#define _nrn_initial _nrn_initial__k_mod_pd
#define nrn_cur _nrn_cur__k_mod_pd
#define _nrn_current _nrn_current__k_mod_pd
#define nrn_jacob _nrn_jacob__k_mod_pd
#define nrn_state _nrn_state__k_mod_pd
#define _net_receive _net_receive__k_mod_pd 
#define rates rates__k_mod_pd 
#define states states__k_mod_pd 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg(int);
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gkbar _p[0]
#define gkbar_columnindex 0
#define A2A_on _p[1]
#define A2A_on_columnindex 1
#define DA_level _p[2]
#define DA_level_columnindex 2
#define pd_cAMP_gain _p[3]
#define pd_cAMP_gain_columnindex 3
#define base_cAMP _p[4]
#define base_cAMP_columnindex 4
#define ik _p[5]
#define ik_columnindex 5
#define g_mod _p[6]
#define g_mod_columnindex 6
#define n_inf _p[7]
#define n_inf_columnindex 7
#define tau_n _p[8]
#define tau_n_columnindex 8
#define n _p[9]
#define n_columnindex 9
#define cAMP _p[10]
#define cAMP_columnindex 10
#define ek _p[11]
#define ek_columnindex 11
#define Dn _p[12]
#define Dn_columnindex 12
#define DcAMP _p[13]
#define DcAMP_columnindex 13
#define v _p[14]
#define v_columnindex 14
#define _g _p[15]
#define _g_columnindex 15
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_vtrap(void);
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
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_k_mod_pd", _hoc_setdata,
 "rates_k_mod_pd", _hoc_rates,
 "vtrap_k_mod_pd", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_k_mod_pd
 extern double vtrap( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define basal_cAMP_prod basal_cAMP_prod_k_mod_pd
 double basal_cAMP_prod = 0.1;
#define deg_rate deg_rate_k_mod_pd
 double deg_rate = 0.1;
#define nq nq_k_mod_pd
 double nq = 4;
#define prod_rate prod_rate_k_mod_pd
 double prod_rate = 0.5;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "prod_rate_k_mod_pd", "/ms",
 "deg_rate_k_mod_pd", "/ms",
 "basal_cAMP_prod_k_mod_pd", "/ms",
 "gkbar_k_mod_pd", "S/cm2",
 "ik_k_mod_pd", "mA/cm2",
 "g_mod_k_mod_pd", "S/cm2",
 "tau_n_k_mod_pd", "ms",
 0,0
};
 static double cAMP0 = 0;
 static double delta_t = 0.01;
 static double n0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "prod_rate_k_mod_pd", &prod_rate_k_mod_pd,
 "deg_rate_k_mod_pd", &deg_rate_k_mod_pd,
 "basal_cAMP_prod_k_mod_pd", &basal_cAMP_prod_k_mod_pd,
 "nq_k_mod_pd", &nq_k_mod_pd,
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
"k_mod_pd",
 "gkbar_k_mod_pd",
 "A2A_on_k_mod_pd",
 "DA_level_k_mod_pd",
 "pd_cAMP_gain_k_mod_pd",
 "base_cAMP_k_mod_pd",
 0,
 "ik_k_mod_pd",
 "g_mod_k_mod_pd",
 "n_inf_k_mod_pd",
 "tau_n_k_mod_pd",
 0,
 "n_k_mod_pd",
 "cAMP_k_mod_pd",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 16, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.036;
 	A2A_on = 0;
 	DA_level = 1;
 	pd_cAMP_gain = 0.5;
 	base_cAMP = 1;
 	_prop->param = _p;
 	_prop->param_size = 16;
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

 void _k_mod_pd_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", 1.0);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 16, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 k_mod_pd k_mod_pd.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   DcAMP = ( A2A_on * ( 1.0 - DA_level ) * prod_rate ) + basal_cAMP_prod - ( deg_rate * cAMP ) ;
   Dn = ( n_inf - n ) / tau_n ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 DcAMP = DcAMP  / (1. - dt*( ( - ( ( deg_rate )*( 1.0 ) ) ) )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_n )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    cAMP = cAMP + (1. - exp(dt*(( - ( ( deg_rate )*( 1.0 ) ) ))))*(- ( ( ( ( A2A_on )*( ( 1.0 - DA_level ) ) )*( prod_rate ) ) + basal_cAMP_prod ) / ( ( - ( ( deg_rate )*( 1.0 ) ) ) ) - cAMP) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_n)))*(- ( ( ( n_inf ) ) / tau_n ) / ( ( ( ( - 1.0 ) ) ) / tau_n ) - n) ;
   }
  return 0;
}
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha_n , _lbeta_n ;
 _lalpha_n = 0.01 * vtrap ( _threadargscomma_ 10.0 - _lv , 10.0 ) ;
   _lbeta_n = 0.125 * exp ( - _lv / 80.0 ) ;
   tau_n = 1.0 / ( _lalpha_n + _lbeta_n ) ;
   n_inf = _lalpha_n / ( _lalpha_n + _lbeta_n ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
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

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  cAMP = cAMP0;
  n = n0;
 {
   rates ( _threadargscomma_ v ) ;
   n = 1.0 ;
   cAMP = base_cAMP ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   rates ( _threadargscomma_ v ) ;
   if ( cAMP < 0.0 ) {
     cAMP = 0.0 ;
     }
   else if ( cAMP > 3.0 ) {
     cAMP = 3.0 ;
     }
   if ( A2A_on  == 1.0 ) {
     g_mod = gkbar * ( 1.0 + pd_cAMP_gain * ( cAMP - base_cAMP ) ) ;
     }
   else {
     g_mod = gkbar ;
     }
   ik = g_mod * pow( n , nq ) * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
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
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = cAMP_columnindex;  _dlist1[0] = DcAMP_columnindex;
 _slist1[1] = n_columnindex;  _dlist1[1] = Dn_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "k_mod_pd.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "This NMODL file defines a potassium current (k_mod_pd) adapted for Parkinson's disease modeling.\n"
  "\n"
  "Enhancements:\n"
  "- Simulates dopamine depletion via DA_level\n"
  "- A2A receptor activation increases cAMP production based on DA level\n"
  "- cAMP dynamically modulates potassium conductance (g_mod)\n"
  "- Optional cAMP kinetics (production/degradation)\n"
  "- MODIFIED: A2A activation now *INCREASES* potassium conductance (facilitatory effect on current)\n"
  "- NEW: Includes a 'basal_cAMP_prod' to maintain a stable baseline cAMP level in control conditions.\n"
  "- TEMPORARY DEBUGGING: 'n' gate is set to 1 to isolate g_mod modulation.\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX k_mod_pd\n"
  "    USEION k READ ek WRITE ik VALENCE 1\n"
  "    RANGE gkbar, ik, g_mod, n_inf, tau_n\n"
  "    RANGE A2A_on, DA_level, base_cAMP, pd_cAMP_gain\n"
  "    GLOBAL nq, prod_rate, deg_rate, basal_cAMP_prod\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (mA) = (milliamp)\n"
  "    (mV) = (millivolt)\n"
  "    (S)  = (siemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    gkbar = 0.036 (S/cm2)       : Max potassium conductance (this is now the BASELINE for modulation)\n"
  "    A2A_on = 0                  : A2A receptor activation toggle (0 = off, 1 = on)\n"
  "    DA_level = 1                : Dopamine level (1 = normal, 0 = PD)\n"
  "    pd_cAMP_gain = 0.5          : Sensitivity of g_mod to cAMP in PD (now represents FACILITATION strength)\n"
  "    base_cAMP = 1               : Baseline cAMP level\n"
  "    prod_rate = 0.5 (/ms)       : A2A-driven cAMP production rate (when A2A_on and DA loss)\n"
  "    deg_rate = 0.1 (/ms)        : cAMP degradation rate\n"
  "    basal_cAMP_prod = 0.1 (/ms) : Basal cAMP production rate (always on)\n"
  "    nq = 4                      : Gating exponent\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    v (mV)\n"
  "    ek (mV)\n"
  "    ik (mA/cm2)\n"
  "    g_mod (S/cm2)\n"
  "    n_inf\n"
  "    tau_n (ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    n                               : Activation gate\n"
  "    cAMP                            : Intracellular modulator (dynamic)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    rates(v)\n"
  "    n = 1.0  : TEMPORARY CHANGE: Set n to always be open for debugging\n"
  "    cAMP = base_cAMP\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE states METHOD cnexp\n"
  "    rates(v)\n"
  "\n"
  "    : Clamp cAMP to safe physiological range\n"
  "    if (cAMP < 0) {\n"
  "        cAMP = 0\n"
  "    } else if (cAMP > 3) {\n"
  "        cAMP = 3\n"
  "    }\n"
  "\n"
  "    : Modulate conductance based on cAMP (with A2A receptor activity)\n"
  "    : A2A activation (via cAMP increase) should *INCREASE* potassium conductance\n"
  "    if (A2A_on == 1) {\n"
  "        : g_mod increases as cAMP increases above base_cAMP\n"
  "        g_mod = gkbar * (1 + pd_cAMP_gain * (cAMP - base_cAMP))\n"
  "        : Optional: Add a check for maximum conductance if needed,\n"
  "        : e.g., if you don't want it to increase indefinitely beyond a certain point.\n"
  "        : Example: if (g_mod > gkbar * 2) { g_mod = gkbar * 2 }\n"
  "    } else {\n"
  "        : In control, A2A_on is 0, so g_mod is just gkbar (no A2A modulation)\n"
  "        g_mod = gkbar\n"
  "    }\n"
  "\n"
  "    ik = g_mod * n^nq * (v - ek)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    : cAMP kinetics: A2A-driven production (when A2A_on and DA loss) + basal production - degradation\n"
  "    cAMP' = (A2A_on * (1 - DA_level) * prod_rate) + basal_cAMP_prod - (deg_rate * cAMP)\n"
  "    n' = (n_inf - n) / tau_n  : Keep n' equation, but n is forced to 1 in INITIAL\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) {\n"
  "    LOCAL alpha_n, beta_n\n"
  "\n"
  "    alpha_n = 0.01 * vtrap(10 - v, 10)\n"
  "    beta_n = 0.125 * exp(-v / 80)\n"
  "\n"
  "    tau_n = 1 / (alpha_n + beta_n)\n"
  "    n_inf = alpha_n / (alpha_n + beta_n)\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x (mV), y (mV)) (1) {\n"
  "    if (fabs(x / y) < 1e-6) {\n"
  "        vtrap = y * (1 - x / y / 2)\n"
  "    } else {\n"
  "        vtrap = x / (exp(x / y) - 1)\n"
  "    }\n"
  "}\n"
  ;
#endif
