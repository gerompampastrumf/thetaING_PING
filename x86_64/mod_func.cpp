#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _CA1ih_reg(void);
extern void _CA1ika_reg(void);
extern void _CA1ikdr_reg(void);
extern void _CA1ina_reg(void);
extern void _ExpGABAab_reg(void);
extern void _Ih_reg(void);
extern void _MyExp2Sid_reg(void);
extern void _MyExp2Sidnw_reg(void);
extern void _MyExp2Syn_reg(void);
extern void _MyExp2SynAlpha_reg(void);
extern void _MyExp2SynBB_reg(void);
extern void _MyExp2SynNMDA_reg(void);
extern void _MyExp2SynNMDABB_reg(void);
extern void _MyExp2SynNMDABB_a_reg(void);
extern void _MyExp2SynNMDABB_b_reg(void);
extern void _SIN_reg(void);
extern void _bgka_reg(void);
extern void _burststim_reg(void);
extern void _caolmw_reg(void);
extern void _ccanl_reg(void);
extern void _ch_CavL_reg(void);
extern void _ch_CavN_reg(void);
extern void _ch_HCN_reg(void);
extern void _ch_HCNolm_reg(void);
extern void _ch_HCNp_reg(void);
extern void _ch_KCaS_reg(void);
extern void _ch_Kdrfast_reg(void);
extern void _ch_Kdrfastngf_reg(void);
extern void _ch_Kdrp_reg(void);
extern void _ch_Kdrslow_reg(void);
extern void _ch_KvA_reg(void);
extern void _ch_KvAdistp_reg(void);
extern void _ch_KvAngf_reg(void);
extern void _ch_KvAolm_reg(void);
extern void _ch_KvAproxp_reg(void);
extern void _ch_KvCaB_reg(void);
extern void _ch_KvGroup_reg(void);
extern void _ch_Nav_reg(void);
extern void _ch_Navaxonp_reg(void);
extern void _ch_Navbis_reg(void);
extern void _ch_Navcck_reg(void);
extern void _ch_Navngf_reg(void);
extern void _ch_Navp_reg(void);
extern void _ch_leak_reg(void);
extern void _constant_reg(void);
extern void _gskch_reg(void);
extern void _icaolmw_reg(void);
extern void _ichan2cck_reg(void);
extern void _iconc_Ca_reg(void);
extern void _iholmw_reg(void);
extern void _kcaolmw_reg(void);
extern void _kdrbwb_reg(void);
extern void _lca_reg(void);
extern void _mykca_reg(void);
extern void _mynetstim_reg(void);
extern void _nafbwb_reg(void);
extern void _nca_reg(void);
extern void _positionfcns_reg(void);
extern void _regn_stim_reg(void);
extern void _vecstim_reg(void);
extern void _xtra_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mechanisms//CA1ih.mod\"");
    fprintf(stderr, " \"mechanisms//CA1ika.mod\"");
    fprintf(stderr, " \"mechanisms//CA1ikdr.mod\"");
    fprintf(stderr, " \"mechanisms//CA1ina.mod\"");
    fprintf(stderr, " \"mechanisms//ExpGABAab.mod\"");
    fprintf(stderr, " \"mechanisms//Ih.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2Sid.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2Sidnw.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2Syn.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynAlpha.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynBB.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynNMDA.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynNMDABB.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynNMDABB_a.mod\"");
    fprintf(stderr, " \"mechanisms//MyExp2SynNMDABB_b.mod\"");
    fprintf(stderr, " \"mechanisms//SIN.mod\"");
    fprintf(stderr, " \"mechanisms//bgka.mod\"");
    fprintf(stderr, " \"mechanisms//burststim.mod\"");
    fprintf(stderr, " \"mechanisms//caolmw.mod\"");
    fprintf(stderr, " \"mechanisms//ccanl.mod\"");
    fprintf(stderr, " \"mechanisms//ch_CavL.mod\"");
    fprintf(stderr, " \"mechanisms//ch_CavN.mod\"");
    fprintf(stderr, " \"mechanisms//ch_HCN.mod\"");
    fprintf(stderr, " \"mechanisms//ch_HCNolm.mod\"");
    fprintf(stderr, " \"mechanisms//ch_HCNp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KCaS.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Kdrfast.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Kdrfastngf.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Kdrp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Kdrslow.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvA.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvAdistp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvAngf.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvAolm.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvAproxp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvCaB.mod\"");
    fprintf(stderr, " \"mechanisms//ch_KvGroup.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Nav.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Navaxonp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Navbis.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Navcck.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Navngf.mod\"");
    fprintf(stderr, " \"mechanisms//ch_Navp.mod\"");
    fprintf(stderr, " \"mechanisms//ch_leak.mod\"");
    fprintf(stderr, " \"mechanisms//constant.mod\"");
    fprintf(stderr, " \"mechanisms//gskch.mod\"");
    fprintf(stderr, " \"mechanisms//icaolmw.mod\"");
    fprintf(stderr, " \"mechanisms//ichan2cck.mod\"");
    fprintf(stderr, " \"mechanisms//iconc_Ca.mod\"");
    fprintf(stderr, " \"mechanisms//iholmw.mod\"");
    fprintf(stderr, " \"mechanisms//kcaolmw.mod\"");
    fprintf(stderr, " \"mechanisms//kdrbwb.mod\"");
    fprintf(stderr, " \"mechanisms//lca.mod\"");
    fprintf(stderr, " \"mechanisms//mykca.mod\"");
    fprintf(stderr, " \"mechanisms//mynetstim.mod\"");
    fprintf(stderr, " \"mechanisms//nafbwb.mod\"");
    fprintf(stderr, " \"mechanisms//nca.mod\"");
    fprintf(stderr, " \"mechanisms//positionfcns.mod\"");
    fprintf(stderr, " \"mechanisms//regn_stim.mod\"");
    fprintf(stderr, " \"mechanisms//vecstim.mod\"");
    fprintf(stderr, " \"mechanisms//xtra.mod\"");
    fprintf(stderr, "\n");
  }
  _CA1ih_reg();
  _CA1ika_reg();
  _CA1ikdr_reg();
  _CA1ina_reg();
  _ExpGABAab_reg();
  _Ih_reg();
  _MyExp2Sid_reg();
  _MyExp2Sidnw_reg();
  _MyExp2Syn_reg();
  _MyExp2SynAlpha_reg();
  _MyExp2SynBB_reg();
  _MyExp2SynNMDA_reg();
  _MyExp2SynNMDABB_reg();
  _MyExp2SynNMDABB_a_reg();
  _MyExp2SynNMDABB_b_reg();
  _SIN_reg();
  _bgka_reg();
  _burststim_reg();
  _caolmw_reg();
  _ccanl_reg();
  _ch_CavL_reg();
  _ch_CavN_reg();
  _ch_HCN_reg();
  _ch_HCNolm_reg();
  _ch_HCNp_reg();
  _ch_KCaS_reg();
  _ch_Kdrfast_reg();
  _ch_Kdrfastngf_reg();
  _ch_Kdrp_reg();
  _ch_Kdrslow_reg();
  _ch_KvA_reg();
  _ch_KvAdistp_reg();
  _ch_KvAngf_reg();
  _ch_KvAolm_reg();
  _ch_KvAproxp_reg();
  _ch_KvCaB_reg();
  _ch_KvGroup_reg();
  _ch_Nav_reg();
  _ch_Navaxonp_reg();
  _ch_Navbis_reg();
  _ch_Navcck_reg();
  _ch_Navngf_reg();
  _ch_Navp_reg();
  _ch_leak_reg();
  _constant_reg();
  _gskch_reg();
  _icaolmw_reg();
  _ichan2cck_reg();
  _iconc_Ca_reg();
  _iholmw_reg();
  _kcaolmw_reg();
  _kdrbwb_reg();
  _lca_reg();
  _mykca_reg();
  _mynetstim_reg();
  _nafbwb_reg();
  _nca_reg();
  _positionfcns_reg();
  _regn_stim_reg();
  _vecstim_reg();
  _xtra_reg();
}

#if defined(__cplusplus)
}
#endif
