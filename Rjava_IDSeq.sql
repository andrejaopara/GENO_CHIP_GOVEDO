SELECT ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL "ID_ZIVALI",
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM zivali ziv
WHERE 
betka.Pasma_Pop_Pv (ziv_id_seq)=1
--ziv.SP1_SIFRA_PASMA in (1)
     and to_char(DAT_ROJSTVO,'yyyy')>2000;

