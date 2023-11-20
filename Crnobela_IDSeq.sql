SELECT ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL "ID_ZIVALI",
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM zivali ziv
;



SELECT ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL "ID_ZIVALI",
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM zivali ziv
WHERE extract ( year from ziv.DAT_ROJSTVO) > 2000;

select * from zivali ziv where  ziv.STEV_ORIG_ZIVAL = '54565740' ;


SELECT ziv.DRZ_ORIG_ZIVAL
  || ziv.STEV_ORIG_ZIVAL "ID_ZIVALI",
  ziv.ZIV_ID_SEQ,
  ziv.DAT_ROJSTVO,
  ziv.SIF_SPOL
FROM zivali ziv
WHERE  ziv.SP1_SIFRA_PASMA=15;

