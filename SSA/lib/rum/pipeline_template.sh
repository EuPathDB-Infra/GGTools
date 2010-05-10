#!/bin/sh
echo "starting..." > OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
BOWTIEEXE -a --best --strata -f GENOMEBOWTIE READSFILE.CHUNK -v 3 --suppress 6,7,8 -p 1 > OUTDIR/X.CHUNK
echo "finished first bowtie run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/make_GU_and_GNU.pl OUTDIR/X.CHUNK OUTDIR/GU.CHUNK OUTDIR/GNU.CHUNK PAIREDEND
echo "finished parsing genome bowtie run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK

BOWTIEEXE -a --best --strata -f TRANSCRIPTOMEBOWTIE READSFILE.CHUNK -v 3 --suppress 6,7,8 -p 1 > OUTDIR/Y.CHUNK
echo "finished second bowtie run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/make_TU_and_TNU.pl OUTDIR/Y.CHUNK GENEANNOTFILE OUTDIR/TU.CHUNK OUTDIR/TNU.CHUNK PAIREDEND
echo "finished parsing transcriptome bowtie run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK

perl SCRIPTSDIR/merge_GU_and_TU.pl OUTDIR/GU.CHUNK OUTDIR/TU.CHUNK OUTDIR/GNU.CHUNK OUTDIR/TNU.CHUNK OUTDIR/BowtieUnique.CHUNK OUTDIR/CNU.CHUNK PAIREDEND
echo "finished merging TU and GU" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/merge_GNU_and_TNU_and_CNU.pl OUTDIR/GNU.CHUNK OUTDIR/TNU.CHUNK OUTDIR/CNU.CHUNK OUTDIR/BowtieNU.CHUNK
echo "finished merging GNU, TNU and CNU" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/make_unmapped_file.pl READSFILE.CHUNK OUTDIR/BowtieUnique.CHUNK OUTDIR/BowtieNU.CHUNK OUTDIR/R.CHUNK PAIREDEND
echo "finished making R" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK

BLATEXE GENOMEBLAT OUTDIR/R.CHUNK -ooc=OOCFILE OUTDIR/R.CHUNK.blat -minScore=MINSCORE -minIdentity=MINIDENTITY SPEED
echo "finished first BLAT run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
MDUSTEXE OUTDIR/R.CHUNK > OUTDIR/R.mdust.CHUNK
echo "finished running mdust on R" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/parse_blat_out.pl OUTDIR/R.CHUNK OUTDIR/R.CHUNK.blat OUTDIR/R.mdust.CHUNK OUTDIR/BlatUnique.CHUNK OUTDIR/BlatNU.CHUNK READLENGTH
echo "finished parsing first BLAT run" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
perl SCRIPTSDIR/merge_Bowtie_and_Blat.pl OUTDIR/BowtieUnique.CHUNK OUTDIR/BlatUnique.CHUNK OUTDIR/BowtieNU.CHUNK OUTDIR/BlatNU.CHUNK OUTDIR/RUM_Unique.CHUNK OUTDIR/RUM_NU.CHUNK PAIREDEND
echo "finished merging Bowtie and Blat" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK

echo "pipeline complete" >> OUTDIR/rum_log.CHUNK
echo `date` >> OUTDIR/rum_log.CHUNK
